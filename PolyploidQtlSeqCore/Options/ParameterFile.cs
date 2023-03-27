
#if DEBUG
using System.Runtime.CompilerServices;

[assembly: InternalsVisibleTo("PolyploidQtlSeqTests")]
[assembly: InternalsVisibleTo("PolyploidQtlSeqWslTests")]
#endif

namespace PolyploidQtlSeqCore.Options
{
    /// <summary>
    /// パラメータファイル
    /// </summary>
    public class ParameterFile
    {
        private const string COMMENT = "#";
        private static readonly char[] _splitter = new[] { '\t' };

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Parameter file.";

        /// <summary>
        /// オプションのShortName
        /// </summary>
        public const string SHORT_NAME = "P";

        /// <summary>
        /// オプションのLongName
        /// </summary>
        public const string LONG_NAME = "paramsFile";

        /// <summary>
        /// パラメータファイルを作成する。
        /// </summary>
        /// <param name="filePath">パラメータファイルのPATH</param>
        public ParameterFile(string filePath)
        {
            Path = filePath;
            IsEnabled = !string.IsNullOrWhiteSpace(filePath) && File.Exists(filePath);
        }

        /// <summary>
        /// パラメータファイルPATHを取得する。
        /// </summary>
        internal string Path { get; }

        /// <summary>
        /// パラメータファイルが有効かどうかを取得する。
        /// </summary>
        internal bool IsEnabled { get; }

        /// <summary>
        /// パラメータファイルを読み込み、LongName => Value辞書を作成する。
        /// </summary>
        /// <param name="toLongNameDictionary">LongName変換辞書</param>
        /// <returns>パラメーター辞書</returns>
        internal IReadOnlyDictionary<string, string> ToParameterDictionary(IReadOnlyDictionary<string, string> toLongNameDictionary)
        {
            var parameterDictionary = new Dictionary<string, string>();
            if (!IsEnabled) return parameterDictionary;

            Console.WriteLine($"Load settings from {Path}.");

            using var reader = new StreamReader(Path);
            while (!reader.EndOfStream)
            {
                var line = reader.ReadLine();
                if (string.IsNullOrEmpty(line)) continue;
                if (line.StartsWith(COMMENT)) continue;

                var items = line.Split(_splitter);
                if (items.Length != 2) continue;

                if (!toLongNameDictionary.TryGetValue(items[0], out var longName))
                    throw new InvalidOperationException($"The {items[0]} option does not exist.");

                if (parameterDictionary.ContainsKey(longName))
                    throw new InvalidOperationException($"{items[0]} is duplicated.");

                parameterDictionary[longName] = items[1];
            }

            return parameterDictionary;
        }
    }
}
