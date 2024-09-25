import os
from pydantic_settings import BaseSettings, SettingsConfigDict
import logging

logging.basicConfig(
    level=logging.DEBUG,  # Set to DEBUG to capture all log levels
    filename="logs.log",
    filemode="w",  # Overwrites the file every time the program runs
    format="%(asctime)s - %(levelname)s - %(message)s"
)


class Settings(BaseSettings):
    DB_HOST: str
    DB_PORT: int
    DB_NAME: str
    DB_USER: str
    DB_PASSWORD: str
    model_config = SettingsConfigDict(
        env_file=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".env")
    )


settings = Settings()


def get_db_url():
    logging.info(f"postgresql+asyncpg://{settings.DB_USER}:{settings.DB_PASSWORD}@"
           f"{settings.DB_HOST}:{settings.DB_PORT}/{settings.DB_NAME}")
    return (
        f"postgresql+asyncpg://{settings.DB_USER}:{settings.DB_PASSWORD}@"
        f"{settings.DB_HOST}:{settings.DB_PORT}/{settings.DB_NAME}"
    )
