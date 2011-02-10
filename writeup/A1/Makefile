TEX=pdflatex
SRC=Readme
TRASH=$(shell ls *.aux *.idx *.log *.toc *.ilg *.ind *.backup)

all:${SRC}
${SRC}:
	${TEX} ${SRC}	
view:all
	evince  ${SRC}.pdf&
clean:
	@for FILE in ${TRASH}; do\
		if [ -e $$FILE ]; then\
			echo "Deleting $$FILE" && rm $$FILE;\
		fi;\
	done

