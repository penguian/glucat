##   -*- autoconf -*-

dnl    This file was part of the KDE libraries/packages
dnl    Copyright (C) 1997 Janos Farkas (chexum@shadow.banki.hu)
dnl              (C) 1997,98,99 Stephan Kulow (coolo@kde.org)
dnl                       modified by Paul Leopardi (paul.leopardi@gmail.com)
dnl                       for GluCat

dnl    This file is free software; you can redistribute it and/or
dnl    modify it under the terms of the GNU Library General Public
dnl    License as published by the Free Software Foundation; either
dnl    version 2 of the License, or (at your option) any later version.

dnl    This library is distributed in the hope that it will be useful,
dnl    but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl    Library General Public License for more details.

dnl    You should have received a copy of the GNU Library General Public License
dnl    along with this library; see the file COPYING.LIB.  If not, write to
dnl    the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
dnl    Boston, MA 02111-1307, USA.

dnl IMPORTANT NOTE:
dnl Please do not modify this file unless you expect your modifications to be
dnl carried into every other module in the repository.

dnl ------------------------------------------------------------------------
dnl Forward compatibility macros (make autoconf 2.13 look like 2.50),
dnl thanks to Raja R Harinath.
dnl ------------------------------------------------------------------------
dnl

# serial 1

ifdef([_AC_PATH_X_XMKMF],[],
   [AC_DEFUN([_AC_PATH_X_XMKMF],[AC_PATH_X_XMKMF])])
ifdef([AC_OUTPUT_SUBDIRS],[],
   [AC_DEFUN([AC_OUTPUT_SUBDIRS],[subdirs=$1; _AC_OUTPUT_SUBDIRS])])

dnl checks originally from acinclude.m4 for KDE:

AC_DEFUN([GLUCAT_CHECK_LIB],
[
  glucat_saved_ldflags="$LDFLAGS"
  LDFLAGS="$LDFLAGS $all_libraries"
  AC_CHECK_LIB($1, $2, $3, $4)
  LDFLAGS="$glucat_saved_ldflags"
])

AC_DEFUN([GLUCAT_CHECK_LIBS],
[
  glucat_saved_ldflags="$LDFLAGS"
  LDFLAGS="$LDFLAGS $all_libraries"
  AC_CHECK_LIB($1, $2, $3, $4, $5)
  LDFLAGS="$glucat_saved_ldflags"
])

AC_DEFUN([GLUCAT_CHECK_HEADER],
[
AC_LANG_SAVE
   glucat_safe_cppflags=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS $all_includes $CXXFLAGS"
   AC_LANG([C++])
   AC_CHECK_HEADER($1, $2, $3, [$4])
   CPPFLAGS=$glucat_safe_cppflags
   AC_LANG_POP([])
])

AC_DEFUN([GLUCAT_CHECK_HEADERS],
[
AC_LANG_SAVE
   glucat_safe_cppflags=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS $all_includes $CXXFLAGS"
   AC_LANG([C++])
   for k_header in $1; do
      AC_CHECK_HEADER($k_header, $2, $3, [$4])
   done
   CPPFLAGS=$glucat_safe_cppflags
   AC_LANG_POP([])
])

AC_DEFUN([GLUCAT_CHECK_CXX11_HEADERS],
[
  save_CXX="$CXX"
  save_CXXFLAGS="$CXXFLAGS"
  save_HAVE_CXX11="$HAVE_CXX11"
  if test -z "$save_HAVE_CXX11"; then
    AX_CXX_COMPILE_STDCXX(11, noext, mandatory)
  fi
  GLUCAT_CHECK_HEADERS($1,
  [
    CXX="$save_CXX"
    CXXFLAGS="$save_CXXFLAGS"
    if test -z "$save_HAVE_CXX11"; then
      AX_CXX_COMPILE_STDCXX(11, noext, mandatory)
    fi
    $2
  ],
  [
    CXX="$save_CXX"
    CXXFLAGS="$save_CXXFLAGS"
    $3
  ],
  [$4])
])

AC_DEFUN([GLUCAT_CHECK_COMPILER_FLAG],
[
AC_MSG_CHECKING([whether $CXX supports -$1])
glucat_cache=`echo $1 | sed 'y% .=/+-%____p_%'`
AC_CACHE_VAL(glucat_cv_prog_cxx_$glucat_cache,
[
AC_LANG_SAVE
  AC_LANG([C++])
  save_CXXFLAGS="$CXXFLAGS"
  CXXFLAGS="$CXXFLAGS -Werror -$1"
  AC_LINK_IFELSE([AC_LANG_PROGRAM([[]], [[ return 0; ]])],[eval "glucat_cv_prog_cxx_$glucat_cache=yes"],[])
  CXXFLAGS="$save_CXXFLAGS"
  AC_LANG_POP([])
])
if eval "test \"`echo '$glucat_cv_prog_cxx_'$glucat_cache`\" = yes"; then
 AC_MSG_RESULT(yes)
 :
 $2
else
 AC_MSG_RESULT(no)
 :
 $3
fi
])

AC_DEFUN([GLUCAT_CHECK_EXTRA_LIBS],
[
AC_MSG_CHECKING([for extra includes])
AC_ARG_WITH(extra-includes, [  --with-extra-includes=DIR
                          adds non standard include paths],
  glucat_use_extra_includes="$withval",
  glucat_use_extra_includes=NONE
)
glucat_extra_includes=
if test -n "$glucat_use_extra_includes" && \
   test "$glucat_use_extra_includes" != "NONE"; then

   ac_save_ifs=$IFS
   IFS=':'
   for dir in $glucat_use_extra_includes; do
     glucat_extra_includes="$glucat_extra_includes $dir"
     USER_INCLUDES="$USER_INCLUDES -I$dir"
   done
   IFS=$ac_save_ifs
   glucat_use_extra_includes="added"
else
   glucat_use_extra_includes="no"
fi
AC_SUBST(USER_INCLUDES)

AC_MSG_RESULT($glucat_use_extra_includes)

glucat_extra_libs=
AC_MSG_CHECKING([for extra libs])
AC_ARG_WITH(extra-libs, [  --with-extra-libs=DIR   adds non standard library paths],
  glucat_use_extra_libs=$withval,
  glucat_use_extra_libs=NONE
)
if test -n "$glucat_use_extra_libs" && \
   test "$glucat_use_extra_libs" != "NONE"; then

   ac_save_ifs=$IFS
   IFS=':'
   for dir in $glucat_use_extra_libs; do
     glucat_extra_libs="$glucat_extra_libs $dir"
     GLUCAT_EXTRA_RPATH="$GLUCAT_EXTRA_RPATH -R $dir"
     USER_LDFLAGS="$USER_LDFLAGS -L$dir"
   done
   IFS=$ac_save_ifs
   glucat_use_extra_libs="added"
else
   glucat_use_extra_libs="no"
fi

AC_SUBST(USER_LDFLAGS)

AC_MSG_RESULT($glucat_use_extra_libs)

])

AC_DEFUN([GLUCAT_CREATE_SUBDIRSLIST],
[

DO_NOT_COMPILE="$DO_NOT_COMPILE admin m4 glucat test_runtime"

if test ! -s $srcdir/subdirs; then
  dnl Note: Makefile.common creates subdirs, so this is just a fallback
  TOPSUBDIRS=""
  files=`cd $srcdir && ls -1`
  dirs=`for i in $files; do if test -d $i; then echo $i; fi; done`
  for i in $dirs; do
    echo $i >> $srcdir/subdirs
  done
fi

if test -s $srcdir/inst-apps; then
  ac_topsubdirs="`cat $srcdir/inst-apps`"
else
  ac_topsubdirs="`cat $srcdir/subdirs`"
fi

for i in $ac_topsubdirs; do
  AC_MSG_CHECKING([if $i should be compiled])
  if test -d $srcdir/$i; then
    install_it="yes"
    for j in $DO_NOT_COMPILE; do
      if test $i = $j; then
        install_it="no"
      fi
    done
  else
    install_it="no"
  fi
  AC_MSG_RESULT($install_it)
  vari=`echo $i | sed -e 's,[[-+]],_,g'`
  if test $install_it = "yes"; then
    TOPSUBDIRS="$TOPSUBDIRS $i"
    eval "$vari""_SUBDIR_included=yes"
  else
    eval "$vari""_SUBDIR_included=no"
  fi
done

AC_SUBST(TOPSUBDIRS)
])
