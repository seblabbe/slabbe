VERSION = 0.1

test: install
	# src/sage/all.py must contain the line "from slabbe import *"
	sage -t --force-lib -i src/slabbe/*

install:
	cd .. && sage -pkg slabbe-$(VERSION)
	cd .. && sage -f slabbe-$(VERSION).spkg
	echo "Effectuer une commande sage pour mettre a jour les path..."
	sage -c "a=randint(1,1000);print 'factor(%s) = %s'%(a,factor(a))"

coverage:
	sage -coverage src/slabbe/*

docbuild:
	sage -docbuild slabbe inventory
	sage -docbuild slabbe html
	sage -docbuild slabbe pdf

