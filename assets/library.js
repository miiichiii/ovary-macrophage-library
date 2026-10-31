'use strict';
(() => {
 const select=document.querySelector('#language');
 const query=new URLSearchParams(location.search).get('lang');
 let saved;try{saved=localStorage.getItem('scrna-offline-language-v1')}catch{}
 function render(lang){
  const language=lang==='en'?'en':'ja';document.documentElement.lang=language;select.value=language;
  document.querySelectorAll('[data-ja][data-en]').forEach(node=>node.textContent=node.dataset[language]);
  document.querySelectorAll('a[data-preserve-language]').forEach(link=>{const url=new URL(link.getAttribute('href'),location.href);url.searchParams.set('lang',language);link.setAttribute('href',url.pathname+url.search+url.hash)});
  try{localStorage.setItem('scrna-offline-language-v1',language)}catch{}
 }
 select.addEventListener('change',()=>render(select.value));render(query==='ja'||query==='en'?query:saved);
})();
