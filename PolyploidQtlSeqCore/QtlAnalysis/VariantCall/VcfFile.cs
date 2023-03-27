namespace PolyploidQtlSeqCore.QtlAnalysis.VariantCall
{
    /// <summary>
    /// VCFファイル
    /// </summary>
    internal class VcfFile
    {
        /// <summary>
        /// VCFファイルを作成する。
        /// </summary>
        /// <param name="filePath">VCFファイルPath</param>
        public VcfFile(string filePath)
        {
            if (string.IsNullOrEmpty(filePath)) throw new ArgumentException(null, nameof(filePath));

            Path = filePath;
        }

        /// <summary>
        /// VCFファイルPathを取得する。
        /// </summary>
        public string Path { get; }

        /// <summary>
        /// Indexファイルを作成する。
        /// </summary>
        /// <returns></returns>
        public async ValueTask CreateIndexFileAsync()
        {
            await Tabix.RunAsync(Path);
        }

        /// <summary>
        /// VCFファイルを圧縮する。
        /// </summary>
        /// <returns>圧縮したVCFファイル</returns>
        public async ValueTask<VcfFile> CompressionAsync()
        {
            return await Bgzip.RunAsync(Path);
        }

        /// <summary>
        /// VCFファイルとIndexファイルを削除する。
        /// </summary>
        public void Delete()
        {
            if (File.Exists(Path)) File.Delete(Path);

            var indexFilePath = Path + ".tbi";
            if (File.Exists(indexFilePath)) File.Delete(indexFilePath);
        }
    }
}
