import logging

import aiogram.types
import aiohttp
from aiogram import Bot, types, Dispatcher
from aiogram.contrib.fsm_storage.memory import MemoryStorage
from aiogram.dispatcher.filters.state import State, StatesGroup
from aiogram.types import InlineKeyboardButton, InlineKeyboardMarkup

from information_networks.config import API_TOKEN, URL

bot = Bot(token=API_TOKEN)
dp = Dispatcher(bot, storage=MemoryStorage())


def _(x): return x


links_data = {}


class ShortenLink(StatesGroup):
    waiting_for_link = State()
    waiting_for_password = State()


@dp.message_handler(commands=['start'])
async def start(message: types.Message):
    await message.answer(
        "Welcome to the URL Shortener bot! Send me a link and, optionally, a password to generate a shortened URL.")


@dp.callback_query_handler(lambda c: c.data.startswith('page'))
async def process_callback(callback_query: types.CallbackQuery):
    user_id = callback_query.from_user.id
    _, page_number, *other = callback_query.data.split(':')

    if not (links := links_data.get(user_id)) or 'update' in other:
        links = await get_links(user_id)
    await send_links_page(links, int(page_number))
    await callback_query.answer()


@dp.message_handler(commands=['links'])
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
    links_per_page = 3
    start_index = (page_number - 1) * links_per_page
    end_index = start_index + links_per_page

    links_slice = links[start_index:end_index]
    links_text = "\n\n".join(
        f"{i}. Clicks: {link['clicks']}\nNew: {link['shortUrl']}\nOld: {link['longUrl']}" for i, link in
        enumerate(links_slice, start_index + 1))

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


# @dp.message_handler(state=ShortenLink.waiting_for_link)
async def process_link(message: types.Message):
    # await ShortenLink.waiting_for_password.set()
    await aiogram.types.ChatActions.typing()
    async with (aiohttp.ClientSession() as session,
                session.post(f"{URL}short-urls",
                             params={"userId": message.chat.id, 'url': message.text}
                             ) as response):
        if response.ok:
            data = await response.json()
        else:
            text = await response.text()
            print(f"Error {text}")
    # Reply with the shortened link
    await message.answer(f"New link:\n{data['shortUrl']}", disable_web_page_preview=True)


async def get_links(user_id):
    async with aiohttp.ClientSession() as session:
        async with session.get(
                f"{URL}short-urls", params={"userId": user_id}) as response:
            if response.ok:
                return await response.json()
            text = await response.text()
            print(f"Error {text}")
    return []


@dp.message_handler(commands=['shorten'])
async def shorten_link(message: types.Message):
    await ShortenLink.waiting_for_link.set()
    await message.reply("Send me the link you want to shorten.")


@dp.message_handler(state="*")
async def incorrect_input(message: types.Message):
    if message.text.startswith('https://'):
        await process_link(message)
    else:
        await message.reply("Please send a valid link.")


if __name__ == '__main__':
    from aiogram import executor

    # from aiogram.contrib.middlewares.logging import LoggingMiddleware
    # dp.middleware.setup(LoggingMiddleware())
    logging.basicConfig(level=logging.INFO)

    executor.start_polling(dp, skip_updates=True)
