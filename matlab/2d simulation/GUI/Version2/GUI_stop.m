global pauseflag record_flag vid

if record_flag
    close(vid);
end

stopflag = true;
pauseflag = false;
axes(hs.ax)
hold off