% % This file will help visualize results by computing the Multidimensional 
% scaling per patient with both Euclidean and Wasserstein distance.


directory = "C:\Users\fanga\OneDrive\Documents\Fall 2023\elec 594"; 
k=1;
j =1;
for i=0:7
    % Import information about Euclidean and Wasserstein distances.
    eucPath = strcat(directory, "\", num2str(i),"Euc.csv");
    wassPath = strcat(directory, "\", num2str(i),"Wass.csv");
    label = strcat(directory, "\", num2str(i),"Labels.csv");

    euc = readmatrix(eucPath);
    wass = readmatrix(wassPath);
    label = readmatrix(label);

    sinus = find(label==0);
    jet = find(label==1);

    [mds1, ~] = cmdscale(euc, 3);
    [mds2, ~] = cmdscale(wass, 3);
    
    if k ==1

        figure;
    end

   
    % Plot the results.
    subplot(4,4, k);
    scatter3(mds1(sinus,1), mds1(sinus,2), mds1(sinus,3), "filled", "DisplayName", 'b');
    hold on
    scatter3(mds1(jet,1), mds1(jet,2), mds1(jet,3), "filled", "DisplayName", 'r');
    hold off
    title("Euclidean MDS Patient = "+num2str(i));
    
    subplot(4,4,k+1);
    scatter3(mds2(sinus,1), mds2(sinus,2), mds2(sinus,3), "filled", "DisplayName", 'b');
    hold on
    scatter3(mds2(jet,1), mds2(jet,2), mds2(jet,3), "filled", "DisplayName", 'r');
    hold off
    title("Wasserstein MDS Patient = "+num2str(i));
    k = k+2;
end