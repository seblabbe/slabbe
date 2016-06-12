all: install test

install:
	# python setup.py install
	#sage -pip install --upgrade -v .
	sage -pip install --upgrade .

develop:
	# python setup.py develop
	sage -pip install --upgrade -e .

old-uninstall:
	rm -fr ~/Applications/sage-git/local/lib/python/site-packages/slabbe*

test: 
	sage -tp --force-lib slabbe/*

coverage:
	sage -coverage src/slabbe/*

docbuild:
	cp -r doc/en/slabbe ~/Applications/sage-git/src/doc/en
	sage -docbuild slabbe inventory
	sage -docbuild slabbe html
	sage -docbuild slabbe pdf
	rm -fr ~/Applications/sage-git/src/sage/doc/en/slabbe*

docs:
	cd docs; sage -sh -c "make html"

register:
	sage -python setup.py sdist
	VERSION=`cat VERSION`; sage -sh -c "twine register dist/slabbe-$$VERSION.tar.gz"

dist:
	sage -python setup.py sdist

upload:
	sage -python setup.py sdist
	VERSION=`cat VERSION`; sage -sh -c "twine upload dist/slabbe-$$VERSION.tar.gz"

publish:
	cp ~/Applications/sage-git/src/doc/output/pdf/en/slabbe/slabbe_ref.pdf ../www/Sage/

.PHONY: all install develop test coverage docbuild docs register dist upload publish
