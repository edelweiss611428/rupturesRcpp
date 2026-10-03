document.addEventListener("DOMContentLoaded", function () {
  // Navbar: the Documentation label links to the Documentation page; its menu opens on
  // hover or keyboard focus instead (see extra.css), so swap pkgdown's toggle button for a link.
  var docsToggle = document.getElementById("dropdown-documentation");
  var docsFirst = docsToggle && docsToggle.parentElement.querySelector(".dropdown-menu a");
  if (docsFirst) {
    var docsLink = document.createElement("a");
    docsLink.id = docsToggle.id;
    docsLink.className = docsToggle.className;
    docsLink.href = docsFirst.getAttribute("href");
    docsLink.setAttribute("aria-haspopup", "true");
    docsLink.textContent = docsToggle.textContent;
    docsToggle.replaceWith(docsLink);
  }

  // Documentation pages: mark the current page in the sidebar and in the navbar, and add
  // previous/next links at the end of the page, in the order the sidebar lists them.
  var chapters = Array.prototype.slice.call(document.querySelectorAll(".guide-nav a.guide-chapter"));
  var main = document.getElementById("main");
  if (!chapters.length || !main) return;

  var navbarLink = document.getElementById("dropdown-documentation");
  if (navbarLink) {
    navbarLink.classList.add("active");
    navbarLink.setAttribute("aria-current", "page");
  }

  var page = function (href) { return href.split("/").pop().replace(/\.html$/, ""); };
  var here = page(window.location.pathname);
  var current = chapters.findIndex(function (a) { return page(a.getAttribute("href")) === here; });
  if (current < 0) return;
  chapters[current].setAttribute("aria-current", "page");

  var link = function (target, cls, label) {
    var a = document.createElement("a");
    a.className = cls;
    a.href = target.getAttribute("href");
    var span = document.createElement("span");
    span.className = "chapter-nav-label";
    span.textContent = label;
    a.appendChild(span);
    a.appendChild(document.createTextNode(target.textContent));
    return a;
  };

  var nav = document.createElement("nav");
  nav.className = "chapter-nav";
  nav.setAttribute("aria-label", "Documentation pages");
  if (current > 0) nav.appendChild(link(chapters[current - 1], "chapter-nav-prev", "Previous"));
  if (current < chapters.length - 1) nav.appendChild(link(chapters[current + 1], "chapter-nav-next", "Next"));
  main.appendChild(nav);
});
