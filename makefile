all: install test

install:
	# python setup.py install
	#sage -pip install --upgrade -v .
	sage -pip install --upgrade .

develop:
	# python setup.py develop
	sage -pip install --upgrade -e .

test: 
	sage -tp --force-lib slabbe/*

coverage:
	sage -coverage slabbe/*

dochtml:
	cd docs; sage -sh -c "make html"

docpdf:
	cd docs; sage -sh -c "make latexpdf"

dist:
	sage -python setup.py sdist

register:
	make dist
	VERSION=`cat VERSION`; sage -sh -c "twine register dist/slabbe-$$VERSION.tar.gz"

upload:
	make dist
	VERSION=`cat VERSION`; sage -sh -c "twine upload dist/slabbe-$$VERSION.tar.gz"

.PHONY: all install develop test coverage dochtml docpdf register dist upload
