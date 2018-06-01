import libtbx.easy_mp

def f(x):
    return x**2

def try_easy_mp():
    result = libtbx.easy_mp.pool_map(func=f,args=range(101))
    print(result)
    return result

class f_caller (object) :
    def __init__ (self) :
        self._obj = "DASH"
    def __call__ (self, x) :
        return try_easy_mp(x, self._obj)