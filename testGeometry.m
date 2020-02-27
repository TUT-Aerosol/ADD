numResults = length(Results);

figure;
for i=1:numResults;
plot(Results{i}.x_ts)
hold on;
xlabel('time')
ylabel('x')
end

figure;
for i=1:numResults;
plot(Results{i}.y_ts)
hold on;
xlabel('time')
ylabel('y')
end

figure;
for i=1:numResults;
plot(Results{i}.z_ts)
hold on;
xlabel('time')
ylabel('z')
end

figure;
for i=1:numResults;
plot(Results{i}.x_ts,Results{i}.y_ts)
hold on;
xlabel('x')
ylabel('y')
end