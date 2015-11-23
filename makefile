VERSION = 0.2

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
	cp /Users/slabbe/Applications/sage-git/src/doc/output/pdf/en/slabbe/slabbe_ref.pdf ../www/Sage/slabbe-$(VERSION).pdf
	cp ../slabbe-$(VERSION).spkg ../www/Sage/
	#git tag -a $(VERSION) -m 'version $(VERSION)' 


