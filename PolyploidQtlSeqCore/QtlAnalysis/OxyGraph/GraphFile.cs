namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ画像ファイル
    /// </summary>
    internal class GraphFile
    {
        /// <summary>
        /// グラフ画像ファイルを作成する。
        /// </summary>
        /// <param name="filePath">画像ファイルのPath</param>
        public GraphFile(string filePath)
        {
            if (string.IsNullOrWhiteSpace(filePath)) throw new ArgumentException(null, nameof(filePath));

            Path = filePath;
        }

        /// <summary>
        /// 画像ファイルのPathをｓ取得する。
        /// </summary>
        public string Path { get; }

        /// <summary>
        /// 画像ファイルを削除する。
        /// </summary>
        public void Delete()
        {
            if (File.Exists(Path)) File.Delete(Path);
        }
    }
}
