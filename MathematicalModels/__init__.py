import logging
default_logger_format = ('%(levelname)s: [%(asctime)s] %(name)s'
                         ' - %(message)s')
default_date_format = '%Y-%m-%d %H:%M:%S'

logging.basicConfig(format=default_logger_format, level=logging.INFO,
                    datefmt=default_date_format)


# __version__ = '1.0.0'
