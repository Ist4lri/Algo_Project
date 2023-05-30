def writing_output_file(filename, model):
    """Writing the output file with command line option

    Args:
        filename (string): the path of the file ot be write.
    """
    with open(filename, "w") as file:

        # Mise en forme des résultats

        # Pour la matrice des scores :
        file.write("Matrice des scores : \n\n")
        for line in model.alignement.output_score_mat:
            file.write("[")
            for element in line:
                file.write(str(element))
                file.write("\t")
            file.write("]\n")
        file.write("\n\n\n")

        # Arbre guide newick :
        file.write("Arbre guide au format newick : \n\n")
        file.write(str(model.alignement.output_newick_upgma))
        file.write("\n\n\n")

        # Alignement multiple :
        file.write("\nAlignement multiple : \n\n")
        for i in range(len(model.alignement.output_clades_names_for_align)):
            file.write(str(model.alignement.output_clades_names_for_align[i]) +
                       " : \t"+str(model.alignement.output_multiple_align[i])+"\n")
        file.write("\n\n\n")

        # Matrice positions conservées :
        file.write("\nMatrice des positions conservées : \n\n")
        for line in model.alignement.output_conserved_pos_mat:
            file.write("[")
            for element in line:
                file.write(str(element))
                file.write("\t")
            file.write("]\n")
        file.write("\n\n\n")

        # Final newick :
        file.write("Arbre newick arpès NJ: \n\n")
        file.write(str(model.alignement.output_newick_final))
        file.write("\n\n\n")
