VERSION = 0.3.b

all: install test

install:
	# python setup.py install
	sage -pip install --upgrade -v .

develop:
	# python setup.py develop
	sage -pip install --upgrade -e .

old-install:
	cd .. && sage -pkg slabbe-$(VERSION)
	cd .. && sage -p slabbe-$(VERSION).spkg
	echo "Effectuer une commande sage pour mettre a jour les path..."
	sage -c "a=randint(1,1000);print 'factor(%s) = %s'%(a,factor(a))"

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

publish:
	cp ~/Applications/sage-git/src/doc/output/pdf/en/slabbe/slabbe_ref.pdf ../www/Sage/slabbe-$(VERSION).pdf
	cp ../slabbe-$(VERSION).spkg ../www/Sage/

# Make my spkg a Python package?
# https://python-packaging.readthedocs.org/en/latest/
# https://python-packaging.readthedocs.org/en/latest/minimal.html
# https://packaging.python.org/en/latest/distributing/
