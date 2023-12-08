import json

import aiogram.types
import aiohttp
from aiogram import Bot, types, Dispatcher
from aiogram.contrib.fsm_storage.memory import MemoryStorage
from aiogram.contrib.middlewares.logging import LoggingMiddleware
from aiogram.dispatcher.filters.state import State, StatesGroup
from aiogram.types import InlineKeyboardButton, InlineKeyboardMarkup

from information_networks.config import API_TOKEN, URL
from information_networks.logger import logger

bot = Bot(token=API_TOKEN)
dp = Dispatcher(bot, storage=MemoryStorage())
dp.middleware.setup(LoggingMiddleware(logger))


def _(x): return x


links_data = {}


class ShortenLink(StatesGroup):
    waiting_for_link = State()
    waiting_for_password = State()


@dp.message_handler(commands=['start'], state='*')
async def start(message: types.Message):
    state = Dispatcher.get_current().current_state()
    await state.reset_state(with_data=False)
    logger.debug(f'Process /start command from {message.chat.id}')

    extra_text = ''
    async with state.proxy() as data:
        user_id = data.get('user_id', None)
        if (args := message.get_args()) and args != user_id:
            logger.info(f'Set new session {args!r} to {message.chat.id}')
            data['user_id'] = args
            extra_text = f'Set new session to you account\n\n'
    return await message.answer(f"{extra_text}Welcome to the URL Shortener bot!\n"
                                f"Send me a link and, optionally, a password to generate a shortened URL.")


@dp.callback_query_handler(lambda c: c.data.startswith('page'))
async def process_callback(callback_query: types.CallbackQuery):
    user_id = callback_query.from_user.id
    _, page_number, *other = callback_query.data.split(':')

    logger.debug(f"Process page for {user_id=} {page_number=}, {other}")

    if not (links := links_data.get(user_id)) or 'update' in other:
        links = await get_links(user_id)
        links_data[user_id] = links
    await send_links_page(links, int(page_number))
    await callback_query.answer()


@dp.message_handler(commands=['links'], state='*')
async def show_links(message: types.Message):
    user_id = message.from_user.id

    await aiogram.types.ChatActions.typing()
    links = await get_links(message.chat.id)
    links_data[user_id] = links
    if links:
        await send_links_page(links, page_number=1)
    else:
        await message.answer(_("You haven't generated any links yet."))


async def send_links_page(links, page_number):
    links = list(reversed(links))
    links_per_page = 5
    start_index = (page_number - 1) * links_per_page
    end_index = start_index + links_per_page

    links_slice = links[start_index:end_index]
    links_text = "\n\n".join(
        f"{i}. Clicks: {link['clicks']}\n"
        f"New: {link['shortUrl']}\n"
        f"Old: {link['longUrl']}" + \
        (f"\nPassword: {link['password']}" if link['password'] else "")
        for i, link in enumerate(links_slice, start_index + 1)
    )

    keyboard = generate_pagination_keyboard(links, page_number, links_per_page)
    if callback := types.CallbackQuery.get_current():
        await callback.message.edit_text(f"Your generated links ({len(links)} total):\n\n{links_text}",
                                         reply_markup=keyboard, disable_web_page_preview=True)
    elif message := types.Message.get_current():
        await message.answer(f"Your generated links ({len(links)} total):\n\n{links_text}",
                             reply_markup=keyboard, disable_web_page_preview=True)


def generate_pagination_keyboard(links, current_page, per_page=3):
    links_count = len(links)
    total_pages = (links_count + per_page - 1) // per_page

    keyboard = InlineKeyboardMarkup(row_width=2)
    keyboard.row(InlineKeyboardButton(_("🔄 Update"), callback_data=f"page:{current_page}:update"))

    buttons = []
    if current_page > 1:
        buttons.append(InlineKeyboardButton(_("◀️ Previous"), callback_data=f"page:{current_page - 1}"))
    if current_page < total_pages:
        buttons.append(InlineKeyboardButton(_("Next ▶️"), callback_data=f"page:{current_page + 1}"))
    keyboard.add(*buttons)

    return keyboard


async def process_link(message: types.Message):
    state = Dispatcher.get_current().current_state()
    logger.debug(f"Process process_link {message}")
    await ShortenLink.waiting_for_password.set()
    async with state.proxy() as data:
        data['api_params'] = dict(url=message.text, user_id=message.chat.id)
        keyboard = InlineKeyboardMarkup()
        keyboard.row(InlineKeyboardButton(_("Skip"), callback_data=f"skip-password"))
        m = await message.answer(f"Send password or press 'Skip' button", reply_markup=keyboard)
        data['message_id'] = m.message_id
    return m


@dp.callback_query_handler(lambda c: c.data == 'skip-password', state=ShortenLink.waiting_for_password)
async def skip_password(call: types.CallbackQuery):
    state = Dispatcher.get_current().current_state()
    logger.debug(f"Process skip_password {call}")

    await ShortenLink.waiting_for_link.set()
    await call.message.edit_text("Password was skipped")
    await aiogram.types.ChatActions.typing()

    async with state.proxy() as data:
        api_params = data['api_params']
    try:
        ok, result = await create_short_url(**api_params, password=None)
    except Exception as e:
        logger.error(f"Unexpected error: {e!r}")
        ok, result = False, repr(e)
    if not ok:
        return await call.message.answer(result)
    await call.message.answer(f"✅🔗: {result}", disable_web_page_preview=True)


@dp.message_handler(state=ShortenLink.waiting_for_password)
async def process_password(message: types.Message):
    state = Dispatcher.get_current().current_state()
    logger.debug(f"Process process_password {message}")
    await ShortenLink.waiting_for_link.set()
    async with state.proxy() as data:
        api_params = data['api_params']
        message_id = data['message_id']
    await aiogram.Bot.get_current().edit_message_reply_markup(message.chat.id, message_id)

    try:
        ok, result = await create_short_url(**api_params, password=message.text)
    except Exception as e:
        logger.error(f"Unexpected error: {e!r}")
        ok, result = False, repr(e)
    if not ok:
        return await message.answer(result)
    await message.answer(f"✅🔗: {result}", disable_web_page_preview=True)


async def create_short_url(url: str, user_id: str, password: str | None):
    params = {"userId": user_id, 'url': url}
    if password:
        params.update(password=password)
    async with (aiohttp.ClientSession() as session,
                session.post(f"{URL}short-urls", params=params) as response):
        logger.info(f"Processed process_link for {user_id}, {response.status}")
        if not response.ok:
            text = await response.text()
            logger.error(f"Error {json.dumps(text)}")
            return False, text
        data = await response.json()
    return True, data['shortUrl']


async def get_links(user_id) -> list:
    state = Dispatcher.get_current().current_state()
    async with state.proxy() as data:
        session_id = data.get('user_id', None)
    return await _get_links(user_id) + ([] if not session_id else await _get_links(session_id))


async def _get_links(user_id):
    logger.debug(f"Process get_links for {user_id}")

    async with aiohttp.ClientSession() as session:
        async with session.get(
                f"{URL}short-urls", params={"userId": user_id}) as response:
            logger.info(f"Processed get_links for {user_id}, {response.status}")
            if response.ok:
                return await response.json()
            text = await response.text()
            logger.error(f"Error on get_links {response.status}. {json.dumps(text)}")
    return []


@dp.message_handler(commands=['shorten'])
async def shorten_link(message: types.Message):
    await ShortenLink.waiting_for_link.set()
    await message.reply("Send me the link you want to shorten.")


@dp.message_handler(state="*")
async def incorrect_input(message: types.Message):
    if not message.text.startswith('https://'):
        return await message.reply("Please send a valid link.")
    return await process_link(message)


async def on_startup(dispatcher: Dispatcher):
    bot = dispatcher.bot
    bot_info = await bot.get_me()
    logger.info(f"Starting bot: {bot_info}")
    await bot.set_my_commands([
        types.BotCommand("start", "Send welcome text"),
        types.BotCommand("links", "View all links"),
    ])


if __name__ == '__main__':
    from aiogram import executor

    executor.start_polling(dp, on_startup=on_startup, skip_updates=False)
