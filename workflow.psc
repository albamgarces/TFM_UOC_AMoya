Funcion variable_de_retorno <- Nombre ( Argumentos )
	
Fin Funcion
Algoritmo dup_annot
	Leer annotation_file
	Si anotation_file = genbank Entonces
		gbf_parser
	SiNo
		acciones_por_falso
	Fin Si
	
FinAlgoritmo
