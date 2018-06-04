import libtbx.easy_mp

def f(x):
    return x**2

def try_easy_mp():
    result = libtbx.easy_mp.pool_map(func=f,args=range(101))
    print(result)
    return result

def f_1(x, magic_string):
    print(magic_string)
    return x**2

def try_f_caller():
    f_call = f_caller(magic_string="CATS")
    result = libtbx.easy_mp.pool_map(fixed_func=f_call, args=range(101))
    print(result)
    return result

class f_caller (object) :
    def __init__ (self, magic_string) :
        self._obj = magic_string
    def __call__ (self, x) :
        return f_1(x,magic_string = self._obj)