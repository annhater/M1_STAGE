# install devtools
install.packages("devtools")

# load devtools
library(devtools)

# install arcdiagram
install_github('arcdiagram', username='gastonstat')

# load arcdiagram
library(arcdiagram)

# Reas contact map files

get_contacts <- function(contacts_file){
    #contacts_file <- contacts_file 
    #residues_A <- c(49, 50, 51)
    
    pdb <- read.csv(contacts_file)
    coords_A <- pdb$atom[pdb$atom$resno %in% residues_A & pdb$atom$elety == "CA" & pdb$atom$chain == chain1, c("x", "y", "z")]
    
    # Extraire les coordonnées du Calpha du résidu 50 de la chaîne B
    coords_B <- pdb$atom[pdb$atom$resno == residue_B & pdb$atom$elety == "CA" & pdb$atom$chain == chain2, c("x", "y", "z")]
    
    # Convertir en matrices/vecteurs
    p1 <- as.numeric(coords_A[1, ])  # Calpha du résidu 49 (chaîne A)
    p2 <- as.numeric(coords_A[2, ])  # Calpha du résidu 50 (chaîne A)
    p3 <- as.numeric(coords_A[3, ])  # Calpha du résidu 51 (chaîne A)
    p  <- as.numeric(coords_B[1, ])  # Calpha du résidu 50 (chaîne B)
    
    # Calculer le vecteur normal au plan défini par p1, p2 et p3
    v1 <- p2 - p1
    v2 <- p3 - p1
    normal <- c(
        v1[2] * v2[3] - v1[3] * v2[2],
        v1[3] * v2[1] - v1[1] * v2[3],
        v1[1] * v2[2] - v1[2] * v2[1]
    )  # Produit vectoriel
    
    normal <- normal / sqrt(sum(normal^2))  # Normalisation du vecteur normal
    
    # Calculer la distance signée entre le point p et le plan (projection orientée)
    distance_signee <- sum((p - p1) * normal)  
    
    # Résultat
    return(distance_signee)
}