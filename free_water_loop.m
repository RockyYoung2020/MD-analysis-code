diary LogFile_2020_8_23;

t_start_all = cputime;
File_now = "全部输出结果\亲水性分析\颗粒球直径30";
cd(File_now);

for ALL_FILE = [5 10 12]
    disp('%------------------------------------------------%');
    fprintf('Now reading folder -> %.0f.\n',ALL_FILE);
    File_cal = strcat("颗粒球直径30-heat",num2str(ALL_FILE),'%');    
    cd(File_cal);% 打开文件夹
    for FILE = [1 0.8 0.5 0.4 0.2 0.08 0.05]; % 循环读取 dump.1%_sorption、dump.2%_sorption、dump.3%_sorption...0.8 0.5 0.4 0.2 0.08 0.05
    t_loop_start=cputime; 
    fprintf('Now reading File -> %.2f.\n',FILE);
    disp('%------------------------------------------------%');
    %C:\Users\YMY\Desktop\无定型吸附模型\吸附过程dump文件\dump.0.03_sorption
    %                                                                                                                                                                                                                                                                   
    str= strcat ('dump.heat_30_',num2str(ALL_FILE),'_',num2str(FILE));
    file = str;

% file is file path of lammps file 
try
    dump = fopen(file,'r');
catch
    error('Dumpfile not found!');
end

i=1;
while feof(dump) == 0
    id = fgetl(dump);
     if (strncmpi(id,'ITEM: TIMESTEP',numel('ITEM: TIMESTEP')))
            timestep(i) = str2num(fgetl(dump));
    else
     if (strncmpi(id,'ITEM: NUMBER OF ATOMS',numel('ITEM: NUMBER OF ATOMS')))
            Natoms(i) = str2num(fgetl(dump));
     else
      if (strncmpi(id,'ITEM: BOX BOUNDS',numel('ITEM: BOX BOUNDS')))
            x_bound(i,:) = str2num(fgetl(dump));
            y_bound(i,:) = str2num(fgetl(dump));
            z_bound(i,:) = str2num(fgetl(dump));
      else
       if (strcmpi(id(1:11),'ITEM: ATOMS'))
            for j = 1 : 1: Natoms
                atom_data(j,:,i) = str2num(fgetl(dump));
            end
            i=i+1;
       end
      end 
     end
   end
end
disp('...Now readng finished...');
%%----------------------------------------------%%

% ITEM: ATOMS id1 type2 x3 y4 z5 vx6 vy7 vz8 fx9 fy10 fz11 v_temp_all12
% ke13 pe14 stress[*]15-20
% type is: 1,2,3,4;
all_frame = 100; % 全部需要计算的帧
bin_num = 20;   % 将体系划分block
%%----------------------------------------------%%
%xl,yl,zl
%%----------------------------------------------%%
xl = x_bound(1,2)-x_bound(1,1);
yl = y_bound(1,2)-y_bound(1,1);
zl = z_bound(1,2)-z_bound(1,1);
%%----------------------------------------------%%
bin_size   = xl/bin_num;          % 体系划分block的大小
Nlocal     = length(atom_data);   % 全部的原子
water_num  = zeros(bin_num,all_frame);
%------------------------------------------------%  
water_cap = zeros(all_frame,1);
water_free_tag=zeros(Nlocal,all_frame);
cutoff = 10;
%------------------------------------------------%  
% save_open= strcat (int2str(FILE),'.txt');
% [fid,message] = fopen(save_open,"wt");
write_data = zeros(Nlocal,6);  %id type x y z tag --> 6
% loop all the frame 
t=cputime;
for frame = 1:all_frame
    t_start = cputime;
    disp('%------------------------------------------------%');
    fprintf('Now the calculating frame is: %d.\n',frame);
%------------------------------------------------%    
    now_frame = atom_data(:,:,frame);
    ID     = now_frame(:,1);
    TYPE   = now_frame(:,2);
    XYZ    = now_frame(:,3:5);
%------------------------------------------------%  
    check_capture_water=0;
    check_cap3 = 0;
    in_loop_i=0;
    out_mark_i=0;
    check_in=0;
    check_34=0;
%------------------------------------------------%   
% here is to slice all bins
    for block = 1:bin_num
        check_waternum=0;
% here is to calculate density/number in each number 
        for i = 1:Nlocal
%------------------------------------------------%
            if(block==1 && (TYPE(i)==3||TYPE(i)==4) )
                check_34 = check_34+1;               
 %------------------------------------------------%               
                    for j = 1:Nlocal                        
                            dx = XYZ(j,1)-XYZ(i,1);
                            dy = XYZ(j,2)-XYZ(i,2);
                            dz = XYZ(j,3)-XYZ(i,3);
                            PBC(dx,dy,dz,xl,yl,zl);
                            distance = sqrt(dx*dx+dy*dy+dz*dz);
 %------------------------------------------------%
                        if(distance<=cutoff && (TYPE(j)==1 || TYPE(j)==2) )%&&ENERGY(i)<=0 
                           check_in = check_in+1;
                           in_loop_i = i;
                           if(check_in==1)
                              out_mark_i =i;
                           end
                           if(in_loop_i~=out_mark_i )
                              %fprintf('the marked i is %d.\n',i) ;
                              check_cap3 =check_cap3+1;
                              check_capture_water = check_capture_water+1;
                              water_free_tag(i,frame)=1;
                           end
                           out_mark_i = i;                                          
                        end                                        
                    end
%------------------------------------------------%                                                        
            end
%------------------------------------------------%
        end        
        water_num(block,frame)= check_waternum;
    end % finish one bin
        water_cap(frame)=check_34 - check_capture_water;% save 

% %------------------------------------------------%  
%     title_first    = "ITEM: TIMESTEP";
%     title_number   = "ITEM: NUMBER OF ATOMS";
%     title_boundary = "ITEM: BOX BOUNDS pp pp pp";
%     title_all      = "ITEM: ATOMS id type x y z watertag";
%         
%     fprintf(fid,'%s\n',title_first);
%     fprintf(fid,'%d\n',timestep(frame));
%     
%     fprintf(fid,'%s\n',title_number);
%     fprintf(fid,'%d\n',Natoms(1));
%     
%     fprintf(fid,'%s\n',title_boundary);
%     fprintf(fid,'%d %d\n',[x_bound(1,1),x_bound(1,2)]);
%     fprintf(fid,'%d %d\n',[y_bound(1,1),y_bound(1,2)]);
%     fprintf(fid,'%d %d\n',[z_bound(1,1),z_bound(1,2)]);     
%         
%     fprintf(fid,'%s\n',title_all);
%     [r,c]=size(write_data(:,:));            % 得到矩阵的行数和列数
%         for i=1:r
%             for j=1:c
%                 fprintf(fid,'%f\t',write_data(i,j));
%             end
%                 fprintf(fid,'\r\n');
%         end 
% %------------------------------------------------% 
%    fprintf('Now the frame %d is finished.\n',frame);
%    %disp('%------------------------------------------------%');
% %------------------------------------------------% 
% t_frame=cputime-t_start;
% fprintf('Now the used time of current frame is: %8.1f s.\n',t_frame);
% end
%------------------------------------------------%  
%------------------------------------------------% 
% t_all=(cputime-t)/60;
% fprintf('Now the used time is: %8.1f mins.\n',t_all);
% fclose(fid);
% disp("-------------------");
% disp("----ALL DONE!!!----");
% disp("-------------------");
% 
% Save Files in Excel
% % ------------------------------------------------% 
% % ------------------------------------------------% 
% %save('water_capture_number','water_cap','-mat');
% xlswrite('water_capture_number.xls',water_cap);
% %save('water_number_each_bin.txt','water_num');
% xlswrite('water_number_each_bin.xls',water_num);
% %save('waternum_mesh.txt','waternum_mesh');
% xlswrite('waternum_mesh.xls',waternum_mesh);
% % ------------------------------------------------% 
% % ------------------------------------------------% 
%  
% Plot contours of water for each bin and save each pic
% % ------------------------------------------------% 
% % ------------------------------------------------% 
% target = griddata(waternum_mesh(:,1),waternum_mesh(:,2),waternum_mesh(:,3),x_in,z_in,'nearest');
% for ii = 1:3:300   
%     contour(waternum_mesh(10:end-10,:,ii)',5);
%     xticklabels({'0.2','0.4','0.6','0.8','1.0',...
%                  '1.2','1.4','1.6','1.8','2.0'})
%     ylim([10 40]);
%     set(gca, 'fontsize', 16);
%     set(gcs,'LineWidth',10);
%     xlabel('X-axis X/D');ylabel('Y-axis Bin Number');
%     saveas(gcf,[num2str(ii),'.jpg']); 
% end
% % ------------------------------------------------% 
% % ------------------------------------------------% 
% fclose(fid);
end
%----------------------------------------------%
R       = 30;
V       = 200*150*150-R^3*3.14*4/3; % A^3
KB         = 1.38*1e-23;% J/K         J/A^3 -> J/M^3 -> pa P = n/V*KB*T

P_water  = ((water_cap./3)/(V*1e-30)*KB*340)';
HUMIDITY = (P_water./2.733/1e2)';
 
NEW_WRITE = zeros(100,1);
NEW_WRITE(:,1) = P_water;
NEW_WRITE(:,2) = HUMIDITY;
 cd ..;
 savefile= strcat (num2str(FILE),'.xlsx');
 make_dir = 0;
 if (make_dir == 0)
    mkdir(strcat("计算结果30-heat",num2str(ALL_FILE),'%'));
    make_dir = make_dir+1;   
 end
 cd(strcat("计算结果30-heat",num2str(ALL_FILE),'%'));
 xlswrite(savefile,NEW_WRITE);
 cd ..;
 cd(File_cal);
%----------------------------------------------%

    end
 cd ..;
end
t_start_stop = (cputime-t_start_all)/60;
fprintf('Now the used time is: %.1f mins.\n',t_start_stop);
disp("-------------------");
disp("----ALL DONE!!!----");
disp("-------------------");
% disp("-----------------------");
% disp("------ALL DONE!!!------");
% disp("-----------------------");
% disp("----Close Computer!----");
% disp("-----------------------");
system('shutdown -s');