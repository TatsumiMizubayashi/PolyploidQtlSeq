using NaturalSort.Extension;

namespace PolyploidQtlSeqCore.IO
{
    /// <summary>
    /// ディレクトリ中にある指定パターンファイルのパスを列挙する。
    /// </summary>
    internal static class FileEnumerator
    {
        private static readonly NaturalSortComparer _naturalSortComparer = StringComparer.OrdinalIgnoreCase.WithNaturalSort();

        /// <summary>
        /// ディレクトリ中にある指定パターンのファイルを並べ替えて列挙する。
        /// </summary>
        /// <param name="dirPath">調査ディレクトリのPATH</param>
        /// <param name="pattern">ファイルパターン</param>
        /// <returns>並べ替え済みファイルパス</returns>
        public static string[] Enumerate(string dirPath, string pattern)
        {
            if (string.IsNullOrWhiteSpace(dirPath)) throw new ArgumentException(null, nameof(dirPath));

            return Directory.EnumerateFiles(dirPath, pattern)
                .OrderBy(x => x, _naturalSortComparer)
                .ToArray();
        }
    }
}
