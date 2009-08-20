# load data
function processfile(filename)
  a = load(filename)
  list = cell()
  last = 1
  for i = 2:size(a)(1)
    if (a(i,1) == -1)
      list{size(list)(2)+1} = a(last:i-1, 2)
      last = i
    endif
  endfor
  list{size(list)(2)+1} = a(last:size(a)(1), 2)
  fflist = cell()
  for l = 1:size(list)(2)
    fflist{l} = fft(list{l})
  endfor
  afflist = cell()
  for l = 1:size(fflist)(2)
    afflist{l} = abs(fflist{l})
  endfor
endfunction
