import { firefox } from 'playwright';

const appUrl = process.env.ATLAS_APP_URL || 'http://127.0.0.1:3838';
const browser = await firefox.launch({ headless: true });
const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
const consoleErrors = [];

page.on('console', (message) => {
  if (message.type() === 'error') consoleErrors.push(message.text());
});
page.on('pageerror', (error) => consoleErrors.push(error.message));

await page.goto(appUrl, { waitUntil: 'domcontentloaded', timeout: 60000 });
await page.waitForSelector('text=Legume Root Nodule Symbiosis Atlas', { timeout: 60000 });
await page.waitForSelector('text=Choose how you want to explore the atlas', { timeout: 60000 });
await page.getByRole('link', { name: /Medicago truncatula/i }).click();
await page.waitForSelector('text=Gene expression', { timeout: 60000 });

if (consoleErrors.length) {
  throw new Error(`Browser console errors:\n${consoleErrors.join('\n')}`);
}

await browser.close();
console.log(`Firefox smoke test passed for ${appUrl}`);
