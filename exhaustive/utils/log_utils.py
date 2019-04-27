import logging

logger = logging.getLogger(__name__)


def log(in_msg=None, out_msg=None):
    def decorator(fn):
        def wrapper(*args, **kwargs):

            logger.debug("Entering {:s}...".format(fn.__name__))

            if in_msg is not None:
                logger.info("{}".format(in_msg))

            result = fn(*args, **kwargs)

            logger.debug("Finished {:s}.".format(fn.__name__))

            if out_msg is not None:
                logger.info("{}".format(out_msg))

            return result

        return wrapper

    return decorator
