program main_PRIMA_model


    !#####################################################################################################
    !PRIMA model V08: modify the Elevation tolerance to be calculated every 1000 iterations (same as WDPM)
    !#####################################################################################################

    use PRIMA_inputs
    use main_PRIMA

    implicit none


    character(len=100)::trash,foutflow
    real, ALLOCATABLE ::exc_water_evap(:)
    integer:: ndays,get_nlines,t
    real::time_start,time_end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !CA model inputs
    character(len=100):: fnam, dir,rep_fnam,fin_wd
    character(len=200)::ascii_data(6)
    integer::rain_flag,ncol,nrow,ntot,inactive_cells,active_cells,n_rvr_cell,method,ras_out,freq_err
    integer:: j,niter,ini_dep_flag,lin_indj
    real*8:: xllcorner,yllcorner
    real, ALLOCATABLE ::CA_outflow(:) !outflow volume
    integer,allocatable::ind_node(:,:),rvr_cell_ind(:)
    real,allocatable::mann(:)!,vel(:),dt(:)
    integer, allocatable::rnk_ind(:),ind_sur(:,:),ind_no_data(:)
    real*8, allocatable:: dem_lin(:),wl_lin(:),wdepth(:)
    real:: time_step,no_data,cell_size,error_thr,vol_thr,errord_thr

    common /dims/ncol,nrow,ntot

    call cpu_time(time_start)



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                     READ exc_water_evap_data                                 !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    fnam='exc_water_evap_data.inp'
    ndays=get_nlines(fnam)
    ndays=ndays-1
    allocate(exc_water_evap(ndays))
    !read excess water data
    open(1, file='exc_water_evap_data.inp',status='old')
    read(1,*)trash
    do j=1,ndays
        read(1,*)trash,exc_water_evap(j)
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                     READ exc_water_evap_data                                 !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                         READ CA model inputs                                 !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call read_PRIMA_inputs(ascii_data,error_thr,errord_thr,vol_thr,rain_flag,ncol,nrow,ntot,inactive_cells,&
                              & active_cells,method,n_rvr_cell,niter,xllcorner,yllcorner,CA_outflow,ind_node,&
                              & rvr_cell_ind,mann,rnk_ind,ind_sur,ind_no_data,dem_lin,wl_lin,time_step,no_data&
                              & ,cell_size,ndays,rep_fnam,dir,ini_dep_flag,ras_out,freq_err)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                     END READ CA model inputs                                 !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !write(*,*)'dem outlet',dem_lin(rvr_cell_ind)
    !run the PRIMA model for all events/days
    do t=1,ndays
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!! RUNING PRIMA model      !!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !add/subtract water
            if (method .ne. 2)then
                !add water depth directly on top of the dem
                !apply the CA algorithm
                wl_lin=wl_lin+(exc_water_evap(t)/1000.0);
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                !search for wl_lin<dem_lin and keep it to dem level (depth=0) (when evaporating water from potholes)
                !for active cells only
                if(exc_water_evap(t) .lt. 0.0)then
                    do j=1,active_cells !loop through active cells only to make sure that there is no cell with water level less than DEM
                        lin_indj=rnk_ind(j)
                        if(wl_lin(lin_indj)<dem_lin(lin_indj))  wl_lin(lin_indj)=dem_lin(lin_indj)
                    end do
                end if
                wl_lin(ind_no_data)=no_data
            end if

            call Run_PRIMA(ascii_data,rain_flag,inactive_cells,active_cells,t,n_rvr_cell,method,niter,&
                            & error_thr,cell_size,no_data,time_step,errord_thr,vol_thr,CA_outflow(t),ind_node,rvr_cell_ind,&
                            & mann,rnk_ind,ind_sur,ind_no_data,dem_lin,wl_lin,rep_fnam,dir,xllcorner,yllcorner,ras_out,freq_err)
            !write(*,*)'CA_outflow',CA_outflow(t)

    end do


    !write the final water distribution at final time step
    allocate(wdepth(ntot))
    wdepth=wl_lin-dem_lin
    if(ras_out.eq. 0) then
    fin_wd='fin_wdepth.asc'
        call write_ascii(wdepth,no_data,xllcorner,yllcorner,cell_size,ind_no_data,active_cells,rnk_ind,&
                            & ascii_data,fin_wd,dir, inactive_cells)
        deallocate(wdepth)
    end if
    !write output
    !%write outflows
    if(method .ne. 1 .and. rain_flag .ne. 0)then
        foutflow=trim(dir)// '/' // 'PRIMA_outflow_cms.txt'
        open(4,file=TRIM(ADJUSTL(foutflow)),status='unknown')
        !write(*,*)'CA_outflow',CA_outflow
        CA_outflow=CA_outflow/time_step
        write(4,'(*(f16.5,/))') CA_outflow(:)
        close(4)
    elseif(method .ne. 1 .and. rain_flag .eq. 0) then
        foutflow=trim(dir)// '/' // 'PRIMA_outflow_m3.txt'
        open(4,file=TRIM(ADJUSTL(foutflow)),status='unknown')
        !write(*,*)'CA_outflow',CA_outflow
        write(4,'(*(f16.5,/))') CA_outflow(:)
        close(4)
    end if
    !write(*,*)CA_outflow
    call cpu_time(time_end)
    !write(*,*)'Elapsed time', (time_end-time_start)/60.0,'min'
    !open(2,file=TRIM(ADJUSTL(rep_fnam)),status='unknown',access = 'append') !create report empty file
    !write(2,*)'Elapsed time', (time_end-time_start)/60.0,'min'
    !close(2) !report file
    stop


end program main_PRIMA_model




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                              Subroutines                                     !!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_nlines(fnam)
    implicit none
    character(len=100)::fnam,str
    integer::n,get_nlines,io

    !write(*,*)fnam
    open(1, file=fnam,status='old')

    n=0
    do
        READ(1,*,IOSTAT=io)  str
        IF (io > 0) THEN
          WRITE(*,*) 'Check input.  Something was wrong'
          EXIT
       ELSE IF (io < 0) THEN
          !WRITE(*,*)  'end of file ', n
          EXIT
       ELSE
          n = n + 1
       END IF
    end do
    close(1)


    get_nlines=n

return
end function get_nlines
