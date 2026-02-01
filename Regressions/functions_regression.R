GMM_rock_v3 = function(c0,c1,R,M,a0,c4){

    .expr_c1 = c1 * log(R)
    .expr_a0 = a0 * R  
    .valueFp = .expr_c1 + .expr_a0
    .valueFm = c4 * M
    .value = c0 + .valueFm + .valueFp

    .grad <- array(0, c(length(.value), 4L), list(NULL, c("c0","c1","a0","c4")))  

    .grad[, "c0"] <- 1.0

    .grad[, "c1"] <- log(R)

    .grad[, "a0"] <- R

    .grad[, "c4"] <- M
    
    attr(.value, "gradient") <- .grad
    .value

    }

get_column_name = function(Period,Fs_Obs = 'Obs'){

    if (Period == -1){
    

    column_name = 'PGV_cm_sec'

    } else if (Period == 0) {
    
        column_name = 'PGA_g'
    
    } else if (floor(Period) == Period && Period != -1 && Period != 0) {
            
        column_name = paste0('T = ',Period, ".0")
    
    } else {

        column_name = paste0('T = ',Period)
    }

    if (Fs_Obs == 'Fs'){
        column_name = paste0('Fs - ',column_name)
    }

    return(column_name)


}

df_Period_3reg = function(database_sinNaN){


    # Inicializar un vector para almacenar índices a eliminar
    rows_to_remove <- c()

    # Iterar sobre cada valor único de NGAsubEQID
    for (j in unique(database_sinNaN$NGAsubEQID)) {
        # Filtrar eventos por NGAsubEQID
        evento <- subset(database_sinNaN, database_sinNaN$NGAsubEQID == j)

        # Aplicar la condición
        if (nrow(evento) < 3) {
        # Almacenar índices de filas a eliminar
        rows_to_remove <- c(rows_to_remove, which(database_sinNaN$NGAsubEQID == j))
        }
    }

    # Eliminar filas del data frame

    if (length(rows_to_remove) > 0){

        database_sinNaN <- database_sinNaN[-rows_to_remove, ]

    }

        return(database_sinNaN)
}

get_R = function(M,Rclst,hnew = 0){

    if (hnew == 0) {

        h1 = -0.82
        h2 = 0.252
        h = 10^(h1+h2*M)
        Rref = sqrt(1+h^2)
        R = sqrt(h^2+Rclst^2)

    } else {
        print('New value of h')
        R = sqrt(hnew^2+Rclst^2)
    }

    return(list(R = R, Rref = Rref))

}
