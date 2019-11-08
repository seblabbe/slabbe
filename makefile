all: install test

install:
	sage -pip install --upgrade --no-index -v .
install3:
	~/GitBox/sage3/sage -pip install --upgrade --no-index -v .

install-internet:
	sage -pip install --upgrade -v .

develop:
	# python setup.py develop
	sage -pip install --upgrade -e .

test: 
	sage -tp --force-lib --show-skipped slabbe/*.py slabbe/*.pyx --log=logs/test.log
	sage -t demos/*.rst
testpython3: 
	~/GitBox/sage3/sage -tp --force-lib --show-skipped slabbe/*.py slabbe/*.pyx --log=logs/test.log
	~/GitBox/sage3/sage -t demos/*.rst

testlong:
	sage -tp --long --force-lib --show-skipped slabbe/*.py slabbe/*.pyx --log=logs/testlong.log

pretestpython3:
	find slabbe -name '*.py' | xargs -n 1 python3 -m py_compile 

coverage:
	sage -coverage slabbe/*

doc:
	cd docs && sage -sh -c "make html"

doc-pdf:
	cd docs && sage -sh -c "make latexpdf"

dist:
	sage -python setup.py sdist

check: dist
	VERSION=`cat VERSION`; sage -sh -c "twine check dist/slabbe-$$VERSION.tar.gz"

upload: dist
	VERSION=`cat VERSION`; sage -sh -c "twine upload dist/slabbe-$$VERSION.tar.gz"

clean: clean-doc

clean-doc:
	cd docs && sage -sh -c "make clean"

.PHONY: all install develop test coverage clean clean-doc doc doc-pdf dist register upload
