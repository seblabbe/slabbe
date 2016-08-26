all: install test

install:
	# python setup.py install
	#sage -pip install --upgrade -v .
	sage -pip install --upgrade --no-index .

develop:
	# python setup.py develop
	sage -pip install --upgrade -e .

test: 
	sage -tp --force-lib slabbe/*.py slabbe/*.pyx

coverage:
	sage -coverage slabbe/*

doc:
	cd docs && sage -sh -c "make html"

doc-pdf:
	cd docs && sage -sh -c "make latexpdf"

dist:
	sage -python setup.py sdist

register: dist
	VERSION=`cat VERSION`; sage -sh -c "twine register dist/slabbe-$$VERSION.tar.gz"

upload: dist
	VERSION=`cat VERSION`; sage -sh -c "twine upload dist/slabbe-$$VERSION.tar.gz"

clean: clean-doc

clean-doc:
	cd docs && sage -sh -c "make clean"

.PHONY: all install develop test coverage clean clean-doc doc doc-pdf dist register upload
