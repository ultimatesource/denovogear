Name: py-editdist
Summary: CPython module to quickly calculate Levenshtein's edit distance
Version: 0.3
Release: 1
Source0: http://www2.mindrot.org/files/py-editdist/py-editdist-%{version}.tar.gz
License: BSD
Group: Development/Libraries
BuildRoot: %{_tmppath}/%{name}-buildroot
Requires: %{__python}
BuildRequires: python-devel, gcc
Url: http://www.mindrot.org/py-editdist.html

%description
"editdist" is a CPython module that calculates the Levenshtein edit
distance between two strings. 

%prep
%setup

%build
%{__python} setup.py build

%install
%{__python} setup.py install --root=$RPM_BUILD_ROOT --record=INSTALLED_FILES
sed -e 's|/[^/]*$||' INSTALLED_FILES | grep "site-packages/" | \
    sort | uniq | awk '{ print "%attr(755,root,root) %dir " $1}' > INSTALLED_DIRS
cat INSTALLED_FILES INSTALLED_DIRS > INSTALLED_OBJECTS

%clean
rm -rf $RPM_BUILD_ROOT

%files -f INSTALLED_OBJECTS
%defattr(-,root,root)
%doc LICENSE README TODO ChangeLog

%changelog
* Wed Jul 05 2006 Damien Miller <djm@mindrot.org>
- Build RPM
