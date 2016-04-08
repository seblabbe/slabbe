VERSION = 0.3.beta

all: install test

test: 
	sage -tp --force-lib src/slabbe/*

install:
	cd .. && sage -pkg slabbe-$(VERSION)
	cd .. && sage -p slabbe-$(VERSION).spkg
	echo "Effectuer une commande sage pour mettre a jour les path..."
	sage -c "a=randint(1,1000);print 'factor(%s) = %s'%(a,factor(a))"

coverage:
	sage -coverage src/slabbe/*

docbuild:
	sage -docbuild slabbe inventory
	sage -docbuild slabbe html
	sage -docbuild slabbe pdf

publish:
	cp ~/Applications/sage-git/src/doc/output/pdf/en/slabbe/slabbe_ref.pdf ../www/Sage/slabbe-$(VERSION).pdf
	cp ../slabbe-$(VERSION).spkg ../www/Sage/

# Make my spkg a Python package?
# https://python-packaging.readthedocs.org/en/latest/
# https://python-packaging.readthedocs.org/en/latest/minimal.html
# https://packaging.python.org/en/latest/distributing/
