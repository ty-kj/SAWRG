%crop and show blue blocks weights ratio

figure;
wc=diag(W);
wwc=reshape(wc,[33,24]);
imshow(wwc*255);
wc_sum=sum(sum(wwc));%default 1

[x,y]=ginput(2);   %先用的ginput函数获取图片中数字的起始坐标
wwc2 =imcrop(wwc,[x(1),y(1),abs(x(1)-x(2)),abs(y(1)-y(2))]);
wc_block=sum(sum(wwc2));

disp(wc_block/wc_sum);%weights ratio
