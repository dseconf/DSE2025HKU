% Build the matrix of transition probabilities for the state of the art
% requires setup.m

setup;

str = repmat(' ',1,8);
for j=nstates:-1:1
  str=sprintf('%s %5.2f',str, cgrid(j));
end

fprintf('%s\n%s\n',str,repmat('-',1,8+6*nstates));

for i = nstates:-1:1

  
  str = sprintf('%5.2f | ', cgrid(i));

  str = sprintf('%s%s', str, repmat(' ',1,6*(nstates-i)));

  for j = i:-1:1
    str=sprintf('%s %5.3f',str, stp(j,i));
  end
  
  fprintf('%s\n%s\n',str,repmat('-',1,8+6*nstates));

end