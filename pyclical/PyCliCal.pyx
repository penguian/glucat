# -*- coding: utf-8 -*-
import types
cdef extern from "PyCliCal.h":
    ctypedef double scalar_t

    ctypedef struct cpp_string "string":
        char* c_str()
    cpp_string *new_string "new string" (char* str)
    void del_string "delete" (cpp_string *str)

    ctypedef struct c_index_set "IndexSet":
        int eq "operator=="(c_index_set rhs)
        c_index_set invert "operator~"()
    c_index_set *new_index_set_from_index_set "new IndexSet" (c_index_set ist)
    c_index_set *new_index_set_from_string "new IndexSet" (char* str)
    c_index_set *new_index_set_from_int "new IndexSet" (int ndx)
    void del_index_set "delete" (c_index_set *mv)
    c_index_set c_or  "operator|"(c_index_set lhs, c_index_set rhs)
    c_index_set c_and "operator&"(c_index_set lhs, c_index_set rhs)
    c_index_set c_xor "operator^"(c_index_set lhs, c_index_set rhs)
    cpp_string index_set_to_repr(c_index_set ist)
    cpp_string index_set_to_str(c_index_set ist)

    ctypedef struct c_clifford "Clifford":
        int eq "operator=="(c_clifford rhs)
        c_clifford neg "operator-"()
        c_clifford call "operator()"(int grade)
        scalar_t getitem "operator[]"(c_index_set ist)
        c_clifford random "random"(c_index_set ist)
        void write(char* msg)
    c_clifford *new_clifford_from_clifford "new Clifford" (c_clifford val)
    c_clifford *new_clifford_from_string "new Clifford" (char* str)
    c_clifford *new_clifford_from_scalar "new Clifford" (scalar_t scr)
    void del_clifford "delete" (c_clifford *mv)
    c_clifford c_add "operator+"(c_clifford lhs, c_clifford rhs)
    c_clifford c_sub "operator-"(c_clifford lhs, c_clifford rhs)
    c_clifford c_mul "operator*"(c_clifford lhs, c_clifford rhs)
    c_clifford c_inn "operator&"(c_clifford lhs, c_clifford rhs)
    c_clifford c_lcn "operator%"(c_clifford lhs, c_clifford rhs)
    c_clifford c_out "operator^"(c_clifford lhs, c_clifford rhs)
    c_clifford c_div "operator/"(c_clifford lhs, c_clifford rhs)
    c_clifford c_inv "inv"(c_clifford mv)
    c_clifford c_involute "involute"(c_clifford mv)
    c_clifford c_reverse "reverse"(c_clifford mv)
    c_clifford c_conj "conj"(c_clifford mv)
    scalar_t c_quad "quad"(c_clifford mv)
    scalar_t c_scalar "scalar"(c_clifford mv)
    c_clifford c_pure "pure"(c_clifford mv)
    c_clifford c_even "even"(c_clifford mv)
    c_clifford c_odd "odd"(c_clifford mv)
    scalar_t c_abs "abs"(c_clifford mv)
    scalar_t c_norm "norm"(c_clifford mv)
    scalar_t c_real "real"(c_clifford mv)
    scalar_t c_imag "imag"(c_clifford mv)
    c_clifford c_pow "pow"(c_clifford mv,int m)
    c_clifford c_sqrt "sqrt"(c_clifford mv)
    c_clifford c_exp "exp"(c_clifford mv)
    c_clifford c_log "log"(c_clifford mv)
    c_clifford c_cos "cos"(c_clifford mv)
    c_clifford c_acos "acos"(c_clifford mv)
    c_clifford c_cosh "cosh"(c_clifford mv)
    c_clifford c_acosh "acosh"(c_clifford mv)
    c_clifford c_sin "sin"(c_clifford mv)
    c_clifford c_asin "asin"(c_clifford mv)
    c_clifford c_sinh "sinh"(c_clifford mv)
    c_clifford c_asinh "asinh"(c_clifford mv)
    c_clifford c_tan "tan"(c_clifford mv)
    c_clifford c_atan "atan"(c_clifford mv)
    c_clifford c_tanh "tanh"(c_clifford mv)
    c_clifford c_atanh "atanh"(c_clifford mv)
    cpp_string clifford_to_repr(c_clifford mv)
    cpp_string clifford_to_str(c_clifford mv)

cdef class index_set:
    cdef c_index_set *thisptr # hold a C++ instance which we're wrapping
    def __cinit__(self, val=0):
        t=type(val)
        if   t==index_set:
            self.thisptr = new_index_set_from_index_set((<index_set>val).thisptr[0])
        elif t==types.IntType:
            self.thisptr = new_index_set_from_int(val)
        elif t==types.StringType:
            self.thisptr = new_index_set_from_string(val)
        else:
            s=repr(val)
            self.thisptr = new_index_set_from_string(s)
    def __dealloc__(self):
        del_index_set(self.thisptr)
    def __richcmp__(self, other, int op):
        cdef int equal
        equal = (<index_set>index_set(self)).thisptr.eq((<index_set>index_set(other)).thisptr[0])
        if op == 2: # ==
            return bool(equal)
        elif op == 3: # !=
            return not equal
        else:
            return NotImplemented
    def __coerce__(self,other):
        return index_set(self),index_set(other)
    def __invert__(self):
        cdef index_set result=index_set()
        result.thisptr[0]=(<index_set>index_set(self)).thisptr.invert()
        return result
    def __and__(self,other):
        cdef index_set result=index_set()
        result.thisptr[0]=c_and((<index_set>index_set(self)).thisptr[0],(<index_set>index_set(other)).thisptr[0])
        return result
    def __or__(self,other):
        cdef index_set result=index_set()
        result.thisptr[0]=c_or((<index_set>index_set(self)).thisptr[0],(<index_set>index_set(other)).thisptr[0])
        return result
    def __xor__(self,other):
        cdef index_set result=index_set()
        result.thisptr[0]=c_xor((<index_set>index_set(self)).thisptr[0],(<index_set>index_set(other)).thisptr[0])
        return result
    def __repr__(self):
        return index_set_to_repr(self.thisptr[0]).c_str()
    def __str__(self):
        return index_set_to_str(self.thisptr[0]).c_str()

cdef class clifford:
    cdef c_clifford *thisptr # hold a C++ instance which we're wrapping
    def __cinit__(self, val=0):
        t=type(val)
        if   t==clifford:
            self.thisptr = new_clifford_from_clifford((<clifford>val).thisptr[0])
        elif t==types.IntType or \
             t==types.LongType or\
             t==types.FloatType:
            self.thisptr = new_clifford_from_scalar(val)
        elif t==types.StringType:
            self.thisptr = new_clifford_from_string(val)
        else:
            s=repr(val)
            self.thisptr = new_clifford_from_string(s)
    def __dealloc__(self):
        del_clifford(self.thisptr)
    def __richcmp__(self, other, int op):
        cdef int equal
        equal = (<clifford>clifford(self)).thisptr.eq((<clifford>clifford(other)).thisptr[0])
        if op == 2: # ==
            return bool(equal)
        elif op == 3: # !=
            return not equal
        else:
            return NotImplemented
    def __coerce__(self,other):
        return clifford(self),clifford(other)
    def __neg__(self):
        cdef clifford result=clifford()
        result.thisptr[0]=(<clifford>clifford(self)).thisptr.neg()
        return result
    def __pos__(self):
        return clifford(self)
    def __add__(self,other):
        cdef clifford result=clifford()
        result.thisptr[0]=c_add((<clifford>clifford(self)).thisptr[0],(<clifford>clifford(other)).thisptr[0])
        return result
    def __radd__(self,other):
        cdef clifford result=clifford()
        result.thisptr[0]=c_add((<clifford>clifford(self)).thisptr[0],(<clifford>clifford(other)).thisptr[0])
        return result
    def __sub__(self,other):
        cdef clifford result=clifford()
        result.thisptr[0]=c_sub((<clifford>clifford(self)).thisptr[0],(<clifford>clifford(other)).thisptr[0])
        return result
    def __rsub__(self,other):
        cdef clifford result=clifford()
        result.thisptr[0]=c_sub((<clifford>clifford(self)).thisptr[0],(<clifford>clifford(other)).thisptr[0])
        return result
    def __mul__(self,other):
        cdef clifford result=clifford()
        result.thisptr[0]=c_mul((<clifford>clifford(self)).thisptr[0],(<clifford>clifford(other)).thisptr[0])
        return result
    def __and__(self,other):
        cdef clifford result=clifford()
        result.thisptr[0]=c_inn((<clifford>clifford(self)).thisptr[0],(<clifford>clifford(other)).thisptr[0])
        return result
    def __mod__(self,other):
        cdef clifford result=clifford()
        result.thisptr[0]=c_lcn((<clifford>clifford(self)).thisptr[0],(<clifford>clifford(other)).thisptr[0])
        return result
    def __or__(self,other):
        cdef clifford result=clifford()
        result.thisptr[0]=c_out((<clifford>clifford(self)).thisptr[0],(<clifford>clifford(other)).thisptr[0])
        return result
    def __xor__(self,other):
        cdef clifford result=clifford()
        result.thisptr[0]=c_out((<clifford>clifford(self)).thisptr[0],(<clifford>clifford(other)).thisptr[0])
        return result
    def __div__(self,other):
        cdef clifford result=clifford()
        result.thisptr[0]=c_div((<clifford>clifford(self)).thisptr[0],(<clifford>clifford(other)).thisptr[0])
        return result
    def __pow__(self,m,dummy):
        cdef clifford result=clifford()
        result.thisptr[0]=c_pow((<clifford>clifford(self)).thisptr[0],m)
        return result
    def __call__(self,grade):
        cdef clifford result=clifford()
        result.thisptr[0]=(<clifford>clifford(self)).thisptr.call(grade)
        return result
    def __getitem__(self,ist):
        return (<clifford>clifford(self)).thisptr.getitem((<index_set>index_set(ist)).thisptr[0])
    def __repr__(self):
        return clifford_to_repr(self.thisptr[0]).c_str()
    def __str__(self):
        return clifford_to_str(self.thisptr[0]).c_str()

def sqrt(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_sqrt((<clifford>clifford(mv)).thisptr[0])
    return result
def inv(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_inv((<clifford>clifford(mv)).thisptr[0])
    return result
def involute(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_involute((<clifford>clifford(mv)).thisptr[0])
    return result
def reverse(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_reverse((<clifford>clifford(mv)).thisptr[0])
    return result
def conj(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_conj((<clifford>clifford(mv)).thisptr[0])
    return result
def quad(mv):
    return c_quad((<clifford>clifford(mv)).thisptr[0])
def scalar(mv):
    return c_scalar((<clifford>clifford(mv)).thisptr[0])
def pure(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_pure((<clifford>clifford(mv)).thisptr[0])
    return result
def even(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_even((<clifford>clifford(mv)).thisptr[0])
    return result
def odd(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_odd((<clifford>clifford(mv)).thisptr[0])
    return result
def abs(mv):
    return c_abs((<clifford>clifford(mv)).thisptr[0])
def norm(mv):
    return c_norm((<clifford>clifford(mv)).thisptr[0])
def real(mv):
    return c_real((<clifford>clifford(mv)).thisptr[0])
def imag(mv):
    return c_imag((<clifford>clifford(mv)).thisptr[0])
def pow(mv,m):
    cdef clifford result=clifford()
    result.thisptr[0]=c_pow((<clifford>clifford(mv)).thisptr[0],m)
    return result
def exp(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_exp((<clifford>clifford(mv)).thisptr[0])
    return result
def log(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_log((<clifford>clifford(mv)).thisptr[0])
    return result
def cos(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_cos((<clifford>clifford(mv)).thisptr[0])
    return result
def acos(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_acos((<clifford>clifford(mv)).thisptr[0])
    return result
def cosh(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_cosh((<clifford>clifford(mv)).thisptr[0])
    return result
def acosh(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_acosh((<clifford>clifford(mv)).thisptr[0])
    return result
def sin(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_sin((<clifford>clifford(mv)).thisptr[0])
    return result
def asin(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_asin((<clifford>clifford(mv)).thisptr[0])
    return result
def sinh(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_sinh((<clifford>clifford(mv)).thisptr[0])
    return result
def asinh(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_asinh((<clifford>clifford(mv)).thisptr[0])
    return result
def tan(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_tan((<clifford>clifford(mv)).thisptr[0])
    return result
def atan(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_atan((<clifford>clifford(mv)).thisptr[0])
    return result
def tanh(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_tanh((<clifford>clifford(mv)).thisptr[0])
    return result
def atanh(mv):
    cdef clifford result=clifford()
    result.thisptr[0]=c_atanh((<clifford>clifford(mv)).thisptr[0])
    return result
def random_clifford(ist):
    cdef clifford result=clifford()
    result.thisptr[0]=(<clifford>clifford(result)).thisptr.random((<index_set>index_set(ist)).thisptr[0])
    return result
