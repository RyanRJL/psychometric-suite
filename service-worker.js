/* Service worker for Psychometric Assistant.
   Strategy: cache-first with network fallback. On first visit the precache
   list is fetched. Subsequent visits load instantly from cache and update in
   the background. Any other requests (e.g. cache-busted assets) are cached
   on first fetch.

   Bump CACHE_VERSION below whenever you ship changes you want to force-
   refresh — old caches will be deleted on activation. */

const CACHE_VERSION = 'psyassist-0.10.5';

/* Bare URLs (no ?v=...) — the SW also caches versioned variants on demand
   via the fetch handler below, so this list is just for first-paint speed. */
const PRECACHE_URLS = [
  './',
  './index.html',
  './styles.css',
  './design-system.css',
  './data.js',
  './app.js',
  './design-system.js',
  './app-effectsize-page.js',
  './app-viz-page.js',
  './app-profile-page.js',
  './manifest.json',
  './icon.svg',
  './icon-maskable.svg'
];

self.addEventListener('install', event => {
  event.waitUntil(
    caches.open(CACHE_VERSION)
      .then(cache => cache.addAll(PRECACHE_URLS).catch(err => {
        /* If any one URL 404s, don't abort the whole install — log and
           continue. The runtime fetch handler will pick up the missing
           file when it's actually requested. */
        console.warn('[sw] partial precache:', err);
      }))
      .then(() => self.skipWaiting())
  );
});

self.addEventListener('activate', event => {
  event.waitUntil(
    caches.keys()
      .then(keys => Promise.all(
        keys.filter(k => k !== CACHE_VERSION).map(k => caches.delete(k))
      ))
      .then(() => self.clients.claim())
  );
});

self.addEventListener('fetch', event => {
  /* Only GET requests are cacheable. Don't intercept anything else. */
  if (event.request.method !== 'GET') return;

  const url = new URL(event.request.url);
  /* Same-origin only — leave Google Fonts and any other CDN to their own
     HTTP caching. */
  if (url.origin !== location.origin) return;

  /* MATCH ON THE FULL URL, QUERY STRING INCLUDED.

     This used to pass { ignoreSearch: true }, which made the ?v= mechanism
     inert: a request for `app.js?v=20260910d` matched the precached bare
     `app.js` and the stale copy was served, whatever the version string said.
     CLAUDE.md calls ?v= and CACHE_VERSION "two separate mechanisms, both
     needed" — with ignoreSearch only the second one did anything, so shipping
     a change meant bumping CACHE_VERSION and hoping, and a ?v= bump alone
     changed nothing at all.

     With it removed, a new ?v= is a cache miss, is fetched, and is cached
     under its own URL. The cost is one extra fetch per file per version, which
     is exactly what a cache-buster is for. The bare URLs stay in the precache
     list for first paint; the navigation fallback below still asks for
     './index.html' bare, which is how it is stored. */
  event.respondWith(
    caches.match(event.request).then(cached => {
      if (cached) return cached;
      return fetch(event.request).then(response => {
        if (response && response.ok){
          const clone = response.clone();
          caches.open(CACHE_VERSION).then(cache => cache.put(event.request, clone));
        }
        return response;
      }).catch(() => {
        /* Offline and not cached — for navigations, fall back to the
           shell so the SPA still loads. */
        if (event.request.mode === 'navigate'){
          return caches.match('./index.html');
        }
        return new Response('Offline', { status: 503, statusText: 'Offline' });
      });
    })
  );
});
