import os
import dotenv
import logging

dotenv.load_dotenv()


logging.getLogger('urllib3').setLevel(logging.ERROR)
logging.getLogger('requests').setLevel(logging.ERROR)


def get_logger(name: str):
    logger = logging.getLogger(name)
    if os.getenv('DEBUG', 'False').lower() == 'true':
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(logging.Formatter(
            '[%(asctime)s] [%(name)s] [%(levelname)s]: %(message)s'
        ))

        file_handler = logging.FileHandler('db_setup.log', mode='w')
        file_handler.setFormatter(logging.Formatter(
            '[%(asctime)s] [%(name)s] [%(levelname)s]: %(message)s'
        ))

        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

    return logger

