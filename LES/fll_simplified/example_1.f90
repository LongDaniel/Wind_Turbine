program example_1

  use fll_mods_m
  implicit none

  type(dnode), pointer :: pnode, ptmp, p_reshaped, pcurr
  type(func_data_set) :: fpar
  logical :: ok
  real(rdouble), pointer :: seg_1(:)
  integer(lint) :: n_list_1, i, j

  !> create an empty list. The length of name is no longer than 16
  pnode => fll_mkdir('Segments', fpar)
  call fll_cat(pnode, 6, .true., fpar)
  ok = fll_write(pnode, 'Segments_stat_1.txt', 100, 'A', fpar)

  !> demand the first child, but don't assign the value
  ptmp => fll_mk('seg_1', 'D', 1_lint, 2_lint, fpar)
  ok = fll_mv(ptmp, pnode, fpar)
  call fll_cat(pnode, 6, .true., fpar)
  ok = fll_write(pnode, 'Segments_stat_2.txt', 100, 'A', fpar)
 
  !> assign the value for first child
  ptmp%d1(:) = (/1.0, 2.0/)
  call fll_cat(pnode, 6, .true., fpar)
  ok = fll_write(pnode, 'Segments_stat_3.txt', 100, 'A', fpar)

  !> put the second child
  ptmp => fll_mk('seg_2', 'D', 1_lint, 2_lint, fpar)
  ptmp%d1(:) = (/3.0, 4.0/)
  ok = fll_mv(ptmp, pnode, fpar)
  call fll_cat(pnode, 6, .true., fpar)
  ok = fll_write(pnode, 'Segments_stat_4.txt', 100, 'A', fpar)

  !> put the third child, and test that two nodes can have the same name
  ptmp => fll_mk('seg_2', 'D', 1_lint, 2_lint, fpar)
  ptmp%d1(:) = (/5.0, 6.0/)
  ok = fll_mv(ptmp, pnode, fpar)
  call fll_cat(pnode, 6, .true., fpar)
  ok = fll_write(pnode, 'Segments_stat_5.txt', 100, 'A', fpar)

  !> count the number of nodes in first list
  n_list_1 = fll_nnodes(pnode, '*', 'D', 1_lint, .FALSE., fpar)
  print *, 'Number of nodes in first list is ', n_list_1  

  !> Create a new list, it is a 2D reshape of first list.
  print *, '---- Reshape the first list and save the result to second list ----'
  p_reshaped => fll_mkdir('Segments_2', fpar)
  ptmp => fll_mk('seg', 'D', n_list_1, 2_lint, fpar)
  pcurr => pnode%pchild
  do i = 1, n_list_1
    !do j = 1, 2
      ptmp%d2(i, :) = pcurr%d1(:)
      pcurr => pcurr%pnext
    !enddo
  enddo
  ok = fll_mv(ptmp, p_reshaped, fpar)
  call fll_cat(p_reshaped, 6, .true., fpar)
  ok = fll_write(p_reshaped, 'Segments_stat_6.txt', 100, 'A', fpar)

  !> delete the first list
  print *, '------------ Delete the first list --------------'
  call fll_rm(pnode, fpar)
  print *, 'Printing the first list. it is expected to be null.'
  call fll_cat(pnode, 6, .true., fpar)
  print *, 'Printing the second list. it is expected to be valid.'
  call fll_cat(p_reshaped, 6, .true., fpar)

  print *, '------------ Delete the second list --------------'
  call fll_rm(p_reshaped, fpar)
  print *, 'Printing the second list. it is expected to be null.'
  call fll_cat(p_reshaped, 6, .true., fpar)
  
  print *, 'All lists have been deleted.'
  
end program example_1
