#!/bin/bash

for file in *.fsp; do
    filename=$(echo "$file" | cut -d'.' -f 1)
    echo $filename
    
    # Elimina las líneas que contienen el carácter "%" y reorganiza los campos.
    sed '/%/d' "$file" | awk -v nx="$nx" -v nz="$nz" '{print $2, $1, $6, $5}' > "$filename.xyz"
    
    # Verifica si la primera fila tiene la misma cantidad de columnas que las filas subsiguientes.
    first_row_column_count=$(head -n 1 "$filename.xyz" | awk '{print NF}')
    subsequent_row_column_count=$(awk 'NR>1 {print NF; exit}' "$filename.xyz")

    # Si la primera fila no tiene la misma cantidad de columnas, elimínala y genera un nuevo archivo ".xyz".
    if [ "$first_row_column_count" -ne "$subsequent_row_column_count" ]; then
        tail -n +2 "$filename.xyz" > "$filename-new.xyz"
        mv "$filename-new.xyz" "$filename.xyz"
    fi
done



