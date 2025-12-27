#! /bin/sh
#
# cvs.sh
#
# This file contains support code from Makefile.common
# It defines a shell function for each known target
# and then does a case to call the correct function.

call_and_fix_autoconf()
{
  $AUTOCONF || exit 1
  if test -r configure.ac.in ; then
    perl -pi -e "print \"if test \\\"x\\\$with_fast_perl\\\" = \\\"xyes\\\"; then\
    \\n  perl -i.bak \\\$ac_aux_dir/conf.change.pl \\\$CONFIG_STATUS\
    \\\\\\n    || mv \\\$CONFIG_STATUS.bak \\\$CONFIG_STATUS\
    \\n  rm -f \\\$CONFIG_STATUS.bak\\nfi\
    \\n\" if /^\\s*chmod\\s+.*\\+x\\s+.*CONFIG_STATUS/;" configure
  fi
}

strip_makefile()
{
  if test -f $makefile_wo; then :; else
    perl -e '$in=0; while ( <> ) { $in = 1 if ($_=~ m/^if /); print $_ unless ($in); $in = 0 if ($_ =~ m/^endif/); }' < Makefile.am.in > $makefile_wo
  fi
}

check_autotool_versions()
{
AUTOCONF_VERSION=`$AUTOCONF --version | head -1`
case $AUTOCONF_VERSION in
  [Aa]utoconf*2.[6-7]* ) : ;;
  "" )
    echo "*** AUTOCONF NOT FOUND!."
    echo "*** GluCat requires autoconf 2.60+"
    exit 1
    ;;
  * )
    echo "*** YOU'RE USING $AUTOCONF_VERSION."
    echo "*** GluCat requires autoconf 2.60+"
    exit 1
    ;;
esac

AUTOHEADER_VERSION=`$AUTOHEADER --version | head -1`
case $AUTOHEADER_VERSION in
  [Aa]utoheader*2.[6-7]* ) : ;;
  "" )
    echo "*** AUTOHEADER NOT FOUND!."
    echo "*** GluCat requires autoheader 2.60+ (part of autoconf)"
    exit 1
    ;;
  * )
    echo "*** YOU'RE USING $AUTOHEADER_VERSION."
    echo "*** GluCat requires autoheader 2.60+ (part of autoconf)"
    exit 1
    ;;
esac

AUTOMAKE_STRING=`$AUTOMAKE --version | head -1`
case $AUTOMAKE_STRING in
  automake*1.1[4-9]* ) : ;;
  unsermake* ) :
    echo "*** YOU'RE USING UNSERMAKE."
    echo "*** GOOD LUCK!! :)"
    ;;
  * )
    echo "*** YOU'RE USING $AUTOMAKE_STRING."
    echo "*** GluCat requires automake 1.14+"
    exit 1
    ;;
esac
}

cvs()
{
check_autotool_versions

### Produce acinclude.m4
if grep '\$(top_srcdir)/acinclude.m4:' $makefile_am >/dev/null; then
  echo "*** Creating acinclude.m4"
  rm -f acinclude.m4 configure.files

  strip_makefile
  $MAKE -f $makefile_wo top_srcdir=. ./acinclude.m4
fi

### Make new subdirs and configure.ac.
### The make calls could be optimized away here,
### with a little thought.
if test -r configure.ac.in; then
  rm -f subdirs configure.ac
  echo "*** Creating list of subdirectories"
  subdirs
  echo "*** Creating configure.ac"
  configure_files
  strip_makefile
  $MAKE -f $makefile_wo top_srcdir=. ./configure.ac || exit 1
fi

echo "*** Creating aclocal.m4"
$ACLOCAL || exit 1
echo "*** Creating configure"
call_and_fix_autoconf

if egrep "^AC_CONFIG_HEADERS" configure.ac >/dev/null 2>&1; then
  echo "*** Creating config.h template"
  $AUTOHEADER || exit 1
fi

echo "*** Creating Makefile templates"
$AUTOMAKE || exit 1

if egrep "^cvs-local:" $makefile_am >/dev/null; then \
  strip_makefile
  $MAKE -f $makefile_wo cvs-local || exit 1
fi

echo "*** Creating date/time stamp"
touch stamp-h.in

echo "*** Finished"
echo "    Don't forget to run ./configure"
echo "    If you haven't done so in a while, run ./configure --help"
}

dist()
{
check_autotool_versions

###
### First build all of the files necessary to do just "make"
###
if grep '\$(top_srcdir)/acinclude.m4:' $makefile_am >/dev/null; then
  strip_makefile
  $MAKE -f $makefile_wo top_srcdir=. ./acinclude.m4
fi
if test -r configure.ac.in; then
  subdirs
  configure_files
  strip_makefile
  $MAKE -f $makefile_wo top_srcdir=. ./configure.ac
fi
$ACLOCAL
$AUTOHEADER
$AUTOMAKE --foreign --include-deps
call_and_fix_autoconf
touch stamp-h.in
if grep "^cvs-local:" $makefile_am >/dev/null; then
  strip_makefile
  $MAKE -f $makefile_wo cvs-local
fi
if grep "^cvs-dist-local:" $makefile_am >/dev/null; then
  strip_makefile
  $MAKE -f $makefile_wo cvs-dist-local
fi
}

subdir_dist()
{
$ACLOCAL
$AUTOHEADER
$AUTOMAKE --foreign --include-deps
call_and_fix_autoconf
}

configure_ac()
{
rm -f configure.ac configure.ac.new
kde_use_qt_param=
test -f configure.files || { echo "need configure.files for configure.ac"; exit 1; }
cat `egrep -v "configure.ac.bot" < configure.files` > configure.ac.new
echo "GLUCAT_CREATE_SUBDIRSLIST" >> configure.ac.new
if test -f Makefile.am.in; then
  subdirs=`cat subdirs`
  for dir in $subdirs; do
    dir=`echo $dir | sed -e "s,[-+],_,g"`
    echo "AM_CONDITIONAL($dir""_SUBDIR_included, test \"x\$$dir""_SUBDIR_included\" = xyes)" >> configure.ac.new
  done
fi
# echo "AC_OUTPUT( \\" >> configure.ac.new
mfs=`find . -type d -print | fgrep -v "/." | \
     sed -e "s#\./##" -e "/^debian/d" | sort`
for i in $mfs; do
  topleveldir=`echo $i| sed -e "s#/.*##"`
  if test -f $topleveldir/configure.ac; then
	continue
  fi
  if test -f $i/Makefile.am; then :; else
	continue
  fi
  if test -s inst-apps; then
    if grep "\"^$topleveldir\"" inst-apps > /dev/null 2>&1; then
	continue
    fi
  fi
  if test "." = "$i"; then
    echo "AC_CONFIG_FILES([ Makefile ])" >> configure.ac.new
  else
    echo "AC_CONFIG_FILES([ $i/Makefile ])" >> configure.ac.new
  fi
  if test -n "$UNSERMAKE"; then
    echo "AC_CONFIG_FILES([ $i/Makefile.rules ])" >> configure.ac.new
  fi
done
egrep '^dnl AC_OUTPUT\(.*\)' `cat configure.files` | sed -e "s#^.*dnl AC_OUTPUT(\(.*\))#AC_CONFIG_FILES([ \1 ])#" >> configure.ac.new
if test -n "$UNSERMAKE"; then
  echo "AC_CONFIG_FILES([ MakeVars ])" >> configure.ac.new
fi
echo "AC_OUTPUT" >> configure.ac.new
modulename=
if test -f configure.ac.in; then
   if head -2 configure.ac.in | egrep "^#MIN_CONFIG\(.*\)$" > /dev/null; then
      kde_use_qt_param=`cat configure.ac.in | sed -n -e "s/#MIN_CONFIG(\(.*\))/\1/p"`
   fi
   if head -2 configure.ac.in | egrep "^#MIN_CONFIG" > /dev/null; then
      line=`grep "^AM_INIT_AUTOMAKE(" configure.ac.in`
      if test -n "$line"; then
	  modulename=`echo $line | sed -e "s#AM_INIT_AUTOMAKE(\([^,]*\),.*#\1#"`
	  VERSION=`echo $line | sed -e "s#AM_INIT_AUTOMAKE([^,]*, *\([^)]*\)).*#\1#"`
      fi
      sed -e "s#AM_INIT_AUTOMAKE([^@].*#dnl PACKAGE set before#" \
          configure.ac.new > configure.ac && mv configure.ac configure.ac.new
   fi
fi
if test -z "$modulename" || test "$modulename" = "@MODULENAME@"; then
   modulename=`pwd`; modulename=`basename $modulename`
fi
if test -z "$VERSION" || test "$VERSION" = "@VERSION@"; then
     VERSION="\"3.0\""
fi
if test -n "$kde_use_qt_param"; then
      sed -e "s#^dnl KDE_USE_QT#KDE_USE_QT($kde_use_qt_param)#" \
      	configure.ac.new > configure.ac && mv configure.ac configure.ac.new
fi
sed -e "s#@MODULENAME@#$modulename#" configure.ac.new |
	sed -e "s#@VERSION@#$VERSION#" > configure.ac
botfiles=`cat configure.files | egrep "configure.ac.bot"`
test -n "$botfiles" && cat $botfiles >> configure.ac
rm -f configure.ac.new
}

configure_files()
{
admindir=NO
for i in . .. ../.. ../../..; do
  if test -x $i/admin; then admindir=$i/admin; break; fi
done
rm -f configure.files
touch configure.files
if test -f configure.ac.in && head -2 configure.ac.in | grep "^#MIN_CONFIG" > /dev/null; then
	echo $admindir/configure.ac.min >> configure.files
fi
test -f configure.ac.in && echo configure.ac.in >> configure.files
list=`find . -name "configure.ac.in" -o -name "configure.ac.bot" | sort`
for i in $list; do if test -f $i && test `dirname $i` != "." ; then
  echo $i >> configure.files
fi; done
test -f configure.ac.mid && echo configure.ac.mid >> configure.files
test -f configure.ac.bot && echo configure.ac.bot >> configure.files
}

subdirs()
{
dirs=
compilefirst=`sed -ne 's#^COMPILE_FIRST[ ]*=[ ]*##p' $makefile_am | head -1`
compilelast=`sed -ne 's#^COMPILE_LAST[ ]*=[ ]*##p' $makefile_am | head -1`
for i in `ls -1`; do
    if test -f $i/Makefile.am; then
	case " $compilefirst $compilelast " in
	  *" $i "*) ;;
	  *) dirs="$dirs $i"
	esac
    fi
done
rm -f _SUBDIRS
for i in $dirs; do
    echo $i >> ./_SUBDIRS
done
if test -f Makefile.am.in; then
  cp Makefile.am.in Makefile.am
  topsubdirs=
  subdirs=
  if test -n "$compilefirst"; then
     topsubdirs='$(COMPILE_FIRST)'
     subdirs='$(COMPILE_FIRST) '
  fi
  if test -n "$UNSERMAKE"; then
    for i in $dirs; do
       vari=`echo $i | sed -e "s,[-+],_,g"`
       echo "if $vari""_SUBDIR_included" >> Makefile.am
       echo "$vari""_SUBDIR=$i" >> Makefile.am
       echo "endif" >> Makefile.am
       topsubdirs="$topsubdirs \$($vari""_SUBDIR)"
    done
  fi
  subdirs="$subdirs"'$(TOPSUBDIRS)'
  if test -n "$compilelast"; then
     topsubdirs="$topsubdirs "'$(COMPILE_LAST)'
     subdirs="$subdirs "'$(COMPILE_LAST)'
  fi
  if test -n "$UNSERMAKE"; then
    echo "SUBDIRS=$topsubdirs" >> Makefile.am
  else
    echo "SUBDIRS=$subdirs" >> Makefile.am
  fi
fi
if test -r subdirs && diff subdirs _SUBDIRS > /dev/null; then
  rm -f _SUBDIRS
fi
test -r _SUBDIRS && mv _SUBDIRS subdirs || true
}

admindir=`echo "$0" | sed 's%[\\/][^\\/][^\\/]*$%%'`
test "x$admindir" = "x$0" && admindir=.

test "x$MAKE" = x && MAKE=make
makefile_am=Makefile.am
makefile_wo=Makefile.am
if test -f Makefile.am.in; then
  makefile_am=Makefile.am.in
  makefile_wo=Makefile.am.in.wo
fi

# Suck in the AUTOCONF detection code
. $admindir/detect-autoconf.sh

###
### Main
###

arg=`echo $1 | tr '\-.' __`
case $arg in
  cvs | dist | subdir_dist | configure_ac | configure_files | subdirs ) $arg ;;
  * ) echo "Usage: cvs.sh <target>"
      echo "Target can be one of:"
      echo "    cvs dist subdir.dist"
      echo "    configure.ac configure.files"
      echo ""
      echo "Usage: anything but $1"
      exit 1 ;;
esac

if test -f Makefile.am.in.wo; then
  rm Makefile.am.in.wo
fi

exit 0
