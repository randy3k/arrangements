all:

clean:
	find . -name '*.so' -exec rm {} \; &
	find . -name '*.o' -exec rm {} \; &
	rm -rf man &
	rm -rf docs

shlib:
	find 'src/next' -name '*.c' | \
		sed -E 's/(.*).c/\1/' | xargs -L1 basename | \
			xargs -L1 -I {} gcc --shared -fPIC -o inst/jl/{}.so src/next/{}.c

document:
	R --slave -e "rmarkdown::render('README.Rmd')" &
	R --slave -e "devtools::document(roclets=c('rd', 'collate', 'namespace'))"

site: document
	R --slave -e "pkgdown::build_site()"

install: document
	R CMD INSTALL --no-multiarch --with-keep.source --install-tests .

check: document
	R --slave -e "devtools::check()"
