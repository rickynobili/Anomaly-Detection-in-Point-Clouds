function displayobjnotexture (obj)
figure()
patch('vertices', obj.v, 'faces', obj.f.v, ...
'FaceVertexCData', rand(8,1));
shading interp
axis square;
axis equal;
view([-37.5, -30]);
end