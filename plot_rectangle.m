function h = plot_rectangle(xmin, xmax, ymin, ymax, formatstring)
  xs = [xmin, xmax, xmax, xmin, xmin];
  ys = [ymin, ymin, ymax, ymax, ymin];
  plot(xs, ys, formatstring)
end