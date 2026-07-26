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

Name:           glucat
Version:        0
Release:        0
Summary:        Library of C++ template classes for Clifford algebras
License:        LGPL-3.0-or-later
Group:          Development/Libraries/C and C++
URL:            https://github.com/penguian/glucat
Source0:        %{name}-%{version}.tar.xz

BuildRequires:  armadillo-devel
BuildRequires:  autoconf
BuildRequires:  automake
BuildRequires:  boost-devel >= 1.66.0
BuildRequires:  doctest-devel
BuildRequires:  eigen3-devel
%if 0%{?suse_version} < 1600
BuildRequires:  gcc13-c++
%else
BuildRequires:  gcc-c++
%endif
BuildRequires:  libtool
BuildRequires:  make
BuildRequires:  perl
BuildRequires:  pkgconfig
BuildRequires:  python-rpm-macros
BuildRequires:  %{python_module Cython}
BuildRequires:  %{python_module devel}
BuildRequires:  %{python_module numpy}
BuildRequires:  %{python_module pytest}
BuildRequires:  %{python_module setuptools}

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

%python_subpackages

%prep
%autosetup -p1 -n %{name}-%{version}

%build
%if 0%{?suse_version} < 1600
export CC=gcc-13
export CXX=g++-13
%endif
make -f admin/Makefile.common bootstrap
%configure --enable-shared=no
%make_build
%make_build -C pyclical

%check
%python_exec -m pytest pyclical/test_pytest_doctests.py

%install
%make_install

%files -n libglucat-devel
%license COPYING
%doc AUTHORS AUTHORS.md ChangeLog DESIGN.md NEWS README README.md TODO TODO.md
%{_includedir}/glucat/
%{_libdir}/pkgconfig/glucat.pc

%files %{python_files}
%{python_sitearch}/PyClical*
%{python_sitearch}/pyclical*

%changelog
