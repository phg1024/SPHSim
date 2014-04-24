function [n, nframes, X] = loadSPHResult(filename)
fid = fopen(filename);

n = fread(fid, 1, 'uint32')
nframes = fread(fid, 1, 'uint32')
scale = fread(fid, 1, 'single')

n*nframes

X = fread(fid, 3*n*nframes, 'single');
size(X)
X = reshape(X, 3*n, nframes);
end