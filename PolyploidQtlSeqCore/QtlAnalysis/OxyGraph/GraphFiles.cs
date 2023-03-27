using SkiaSharp;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフファイルコレクション
    /// </summary>
    internal class GraphFiles
    {
        private readonly GraphFile[] _graphFiles;

        /// <summary>
        /// グラフファイルコレクションを作成する。
        /// </summary>
        /// <param name="files"></param>
        public GraphFiles(GraphFile[] files)
        {
            _graphFiles = files;
        }

        /// <summary>
        /// グラフ画像を縦に連結した画像を作成する。
        /// </summary>
        /// <param name="mergeImageFilePath">連結画像ファイルのPath</param>
        /// <param name="option">グラフオプション</param>
        public void VerticalMerge(string mergeImageFilePath, GraphOption option)
        {
            if (string.IsNullOrWhiteSpace(mergeImageFilePath)) throw new ArgumentException(null, nameof(mergeImageFilePath));

            var mergeImageWidth = option.FigureWidth.Value;
            var mergeImageHeight = option.FigureHeight.Value * _graphFiles.Length;
            using var mergeBitmap = new SKBitmap(mergeImageWidth, mergeImageHeight);
            using var canvas = new SKCanvas(mergeBitmap);
            var y = 0;

            foreach (var graphFile in _graphFiles)
            {
                using var graphBitmap = SKBitmap.Decode(graphFile.Path);
                canvas.DrawBitmap(graphBitmap, 0, y);
                y += graphBitmap.Height;
            }

            using var pngData = mergeBitmap.Encode(SKEncodedImageFormat.Png, 100);
            using var stream = File.Open(mergeImageFilePath, FileMode.Create, FileAccess.Write);
            pngData.SaveTo(stream);
        }

        /// <summary>
        /// グラフ画像を削除する。
        /// </summary>
        public void Delete()
        {
            foreach (var graphFile in _graphFiles)
            {
                graphFile.Delete();
            }
        }
    }
}
