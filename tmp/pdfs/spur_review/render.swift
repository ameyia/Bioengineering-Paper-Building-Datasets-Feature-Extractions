import Foundation
import PDFKit
import AppKit

let args = CommandLine.arguments
guard args.count == 3 else { fatalError("usage: render input.pdf output-dir") }
let input = URL(fileURLWithPath: args[1])
let out = URL(fileURLWithPath: args[2], isDirectory: true)
try FileManager.default.createDirectory(at: out, withIntermediateDirectories: true)
guard let doc = PDFDocument(url: input) else { fatalError("cannot open PDF") }
print("PAGES \(doc.pageCount)")
for index in 0..<doc.pageCount {
    guard let page = doc.page(at: index) else { continue }
    print("\n=== SLIDE \(index + 1) ===\n\(page.string ?? "")")
    let box = page.bounds(for: .mediaBox)
    let scale: CGFloat = 1.7
    let size = NSSize(width: box.width * scale, height: box.height * scale)
    let image = NSImage(size: size)
    image.lockFocus()
    NSColor.white.setFill()
    NSRect(origin: .zero, size: size).fill()
    guard let context = NSGraphicsContext.current?.cgContext else { fatalError("no context") }
    context.saveGState()
    context.scaleBy(x: scale, y: scale)
    page.draw(with: .mediaBox, to: context)
    context.restoreGState()
    image.unlockFocus()
    guard let tiff = image.tiffRepresentation,
          let rep = NSBitmapImageRep(data: tiff),
          let png = rep.representation(using: .png, properties: [:]) else { continue }
    try png.write(to: out.appendingPathComponent(String(format: "slide-%02d.png", index + 1)))
}
