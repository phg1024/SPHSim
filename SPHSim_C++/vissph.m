clear all;

[n, nframes, X] = loadSPHResult('logfile.bin');

for i=1:n
    pos = X(:, i);
    pos = reshape(pos, 3, n);
    plot(pos(1,:), pos(2,:), 'b.');
    axis([0 1 0 1]);
    M(i) = getframe;
end

movie(M, 30);