module disulfide_module
    implicit none
    contains

    subroutine get_disulfide_bonds_matrix(pdb_code, SS_matrix, num_rows, num_cols)
        implicit none
        character(len=*), intent(in) :: pdb_code
        integer, allocatable, intent(out) :: SS_matrix(:,:)
        integer, intent(out) :: num_rows, num_cols
        integer :: i, j, ios
        character(len=100) :: line

        ! Llamar al script Python con el código PDB
        call system("python SSbridge.py " // trim(pdb_code), ios)

        ! Verificar si el script Python se ejecutó correctamente
        if (ios /= 0) then
            print *, "Error: Python script execution failed."
            stop
        endif

        ! Abrir el archivo CSV para contar las filas (pares de disulfuro)
        open(unit=10, file='disulfide_bonds.csv', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening the CSV file."
            stop
        endif

        ! Inicializar el contador de filas
        num_rows = 0
        num_cols = 2  ! Sabemos que cada fila tiene 2 columnas (pares de residuos)

        ! Contar las líneas del archivo (número de filas en el archivo CSV)
        do
            read(10, '(A)', iostat=ios) line
            if (ios /= 0) exit
            num_rows = num_rows + 1
        end do

        ! Cerrar y reabrir el archivo para leer los datos después de contar las filas
        close(10)
        open(unit=10, file='disulfide_bonds.csv', status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error reopening the CSV file."
            stop
        endif

        ! Asignar memoria para la matriz dinámica
        allocate(SS_matrix(num_rows, num_cols))

        ! Leer los valores desde el archivo CSV y almacenarlos en la matriz
        do i = 1, num_rows
            read(10, *) (SS_matrix(i, j), j=1,num_cols)
        end do

        close(10)

    end subroutine get_disulfide_bonds_matrix

end module disulfide_module
