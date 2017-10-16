;extend.pro
;by Ken Liao, based on CMREPLICATE, by Craig Markwardt
;
; PURPOSE:
;   Extend an array, increasing its dimensions, duplicating entries across the
;   new dimensions. Behavior is similar to cmreplicate, but extends the
;   dimension to the LEFT instead of the RIGHT.
;
; SYNTAX:
;   RARRAY = EXTEND(DIMS, ARRAY)
;
; example:
;A is a size [2,3] matrix
;extend([5,4],A) is a size [5,4,2,3] array
;
; BUGS:
;   Does not work correctly for input arrays of non-numeric types.
;   But the main point is to use with numeric types, so maybe ok.
;
;CMREPLICATE
; Copyright (C) 2000, 2009, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

function extend, dims, array

  if n_params() EQ 0 then begin
      message, 'RARRAY = EXTEND(DIMS, ARRAY)', /info
      return, 0L
  endif
  on_error, 2

  if n_elements(dims) EQ 0 then return, array
  if n_elements(array) EQ 0 then $
    message, 'ERROR: ARRAY must have at least one element'

  ;; Construct new dimensions, being careful about scalars
  sz = size(array)
  type = sz(sz(0)+1)
  if sz(0) EQ 0 then return, make_array(value=array, dimension=dims)
  onedims = [dims*0+1, sz(1:sz(0))] ;; For REFORM, to extend # of dims.
  newdims = [dims,     sz(1:sz(0))] ;; For REBIN, to enlarge # of dims.
  nnewdims = n_elements(newdims)
  
  if nnewdims GT 8 then $
    message, 'ERROR: resulting array would have too many dimensions.'

  if type NE 7 AND type NE 8 AND type NE 10 AND type NE 11 then begin
      ;; Handle numeric types

      ;; Passing REBIN(...,ARRAY) is a feature introduced in
      ;; IDL 5.6, and will be incompatible with previous versions.
      array1 = rebin(reform([array], onedims), newdims)

      return, array1

  endif else begin
      ;; Handle strings, structures, pointers, and objects separately

      ;; Handle structures, which are never scalars
      if type EQ 8 AND sz(0) EQ 1 AND n_elements(array) EQ 1 then $
        return, reform(make_array(value=array, dimension=dims), dims, /over)

      nold = n_elements(array)
      nadd = 1L
      for i = 0L, n_elements(dims)-1 do nadd = nadd * round(dims(i))
      if type EQ 8 then $
        array1 = make_array(value=array(0), dimension=[nadd,nold]) $
      else $
        array1 = make_array(type=type, dimension=[nadd,nold], /nozero)
      array1 = reform(array1, [nadd,nold], /overwrite)
      array2 = reform([array], n_elements(array))

      ;; Efficient copying, done by powers of two
      array1(0,0) = temporary(array2)
      stride = 1L   ;; stride increase by a factor of two in each iteration
      i = 1L & nleft = nadd - 1
      while nleft GT stride do begin
          array1(i,0) = array1(0:stride-1,*)  ;; Note sneaky IDL optimization
          i = i + stride & nleft = nleft - stride
          stride = stride * 2
      endwhile
      if nleft GT 0 then array1(i,0) = array1(0:nleft-1,*)

      return, reform(array1, newdims, /overwrite)
  endelse

end


