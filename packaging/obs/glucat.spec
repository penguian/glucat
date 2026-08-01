#
# spec file for package glucat
#
# Copyright (c) 2026 SUSE LLC
# Copyright (c) 2026 Paul Leopardi <paul.leopardi@anu.edu.au>
#
# All modifications and additions to the file contributed by third parties
# remain the property of their respective owners, unless otherwise agreed
# upon. The license to this file is the same license as for the pristine
# package itself (Universal Permissive License v 1.0 or LGPL-3.0-or-later).
#

%global flavor @BUILD_FLAVOR@%{nil}

%define pname glucat

%if "%{flavor}" == "main"
%bcond_with    python
%define psuffix -main
%else
%bcond_without python
%define skip_python2 1
%endif

Name:           %{pname}%{?psuffix}
Version:        0
Release:        0
Summary:        Library of C++ template classes for Clifford algebras
License:        LGPL-3.0-or-later
Group:          Development/Libraries/C and C++
URL:            https://github.com/penguian/glucat
Source0:        %{pname}-%{version}.tar.xz

BuildRequires:  autoconf
BuildRequires:  automake
BuildRequires:  boost-devel >= 1.68.0
BuildRequires:  eigen3-devel
BuildRequires:  libtool
BuildRequires:  make
BuildRequires:  perl
BuildRequires:  pkgconfig
%if 0%{?suse_version} < 1600
BuildRequires:  gcc14-c++
%else
BuildRequires:  gcc-c++
%endif

%if %{with python}
BuildRequires:  python-rpm-macros
BuildRequires:  %{python_module Cython}
BuildRequires:  %{python_module devel}
BuildRequires:  %{python_module numpy}
BuildRequires:  %{python_module pytest}
BuildRequires:  %{python_module setuptools}
%python_subpackages
%else
BuildRequires:  armadillo-devel
BuildRequires:  doctest-devel
%endif

%description
GluCat is a library of C++ template classes for Clifford algebras
over the field of real numbers.

%package -n libglucat-devel
Summary:        Header-only C++ template library for Clifford algebras
Group:          Development/Libraries/C and C++
Provides:       glucat-devel = %{version}

%description -n libglucat-devel
GluCat is a library of C++ template classes for Clifford algebras
over the field of real numbers.

This package contains the C++ header files and documentation.

%prep
%autosetup -p1 -n %{pname}-%{version}

%build
%if 0%{?suse_version} < 1600
export CC=gcc-14
export CXX=g++-14
%endif
make -f admin/Makefile.common bootstrap
# Strip -march=native for reproducible OBS builds
sed -i "s|-march=native||g" configure

%if %{with python}
%{python_expand # Apply to all supported python flavors
export PYTHON=$python
mkdir ../${PYTHON}_build
cp -pr ./ ../${PYTHON}_build
pushd ../${PYTHON}_build
%configure --enable-shared=no --enable-pyclical
%make_build
popd
}
%else
%configure --enable-shared=no --disable-pyclical
%make_build
%endif

%check
%if %{with python}
%{python_expand # Apply to all supported python flavors
export PYTHON=$python
pushd ../${PYTHON}_build
$python -m pytest pyclical/test_pytest_doctests.py
popd
}
%endif

%install
%if %{with python}
%{python_expand # Apply to all supported python flavors
export PYTHON=$python
pushd ../${PYTHON}_build
%make_install
popd
}
# Only the compiled Python extension belongs in this flavor.
# Docs and demos are installed by the main flavor or deferred to glucat-doc.
rm -rf %{buildroot}%{_includedir}/
rm -rf %{buildroot}%{_datadir}/doc/
rm -rf %{buildroot}%{_docdir}/
rm -rf %{buildroot}%{_datadir}/pyclical/
%else
%make_install
rm -rf %{buildroot}%{_datadir}/doc/glucat
%endif

%if %{without python}
%files -n libglucat-devel
%license COPYING
%doc AUTHORS AUTHORS.md ChangeLog DESIGN.md NEWS README README.md TODO TODO.md
%{_includedir}/glucat/
%endif

%if %{with python}
%files %{python_files}
%{python_sitearch}/PyClical*.so
%endif

%changelog
