// script.js
window.addEventListener('load', function() {
  var images = document.querySelectorAll('.casID_sprite');
  var maxHeight = 200;

  images.forEach(function(image) {
    var aspectRatio = image.naturalWidth / image.naturalHeight;
    var newWidth = maxHeight * aspectRatio;
    image.style.maxHeight = maxHeight + 'px';
    image.style.width = newWidth + 'px';
  });
});
