def append_docstring(func):
    def decorator(self):
        if self.__doc__:
            self.__doc__ += "\n\n" + func
        else:
            self.__doc__ = func
        return self
    return decorator