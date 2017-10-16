;print time since tic was called
pro toc
  common tictoc_var, t

  t2 = systime(1) - t
  print, t2, format="('time elapsed: ',g0.3,' s')"
end
