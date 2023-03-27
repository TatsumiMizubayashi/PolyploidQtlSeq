namespace PolyploidQtlSeqCore.IO
{
    /// <summary>
    /// BAMファイル列挙
    /// </summary>
    internal static class BamFileEnumerator
    {
        private const string BAM_FILE_PATTERN = "*.bam";

        /// <summary>
        /// BAMファイル Pathを列挙する。
        /// </summary>
        /// <param name="dirPath">調査ディレクトリPath</param>
        /// <returns>BAMファイルPath配列</returns>
        public static string[] Enumerate(string dirPath)
        {
            return FileEnumerator.Enumerate(dirPath, BAM_FILE_PATTERN);
        }
    }
}
