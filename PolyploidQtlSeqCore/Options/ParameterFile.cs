namespace PolyploidQtlSeqCore.Options
{
    /// <summary>
    /// パラメーターファイル
    /// </summary>
    public sealed class ParameterFile
    {
        /// <summary>
        /// オプションのShortName
        /// </summary>
        public const string SHORT_NAME = "P";

        /// <summary>
        /// オプションのLongName
        /// </summary>
        public const string LONG_NAME = "paramsFile";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Parameter file.";

        /// <summary>
        /// コメントアウトマーク
        /// </summary>
        private const char _commentOut = '#';

        /// <summary>
        /// デリミタ
        /// </summary>
        private static readonly char _delimiter = '\t';

        /// <summary>
        /// Key index
        /// </summary>
        private static readonly int _keyIndex = 0;

        /// <summary>
        /// Value index
        /// </summary>
        private static readonly int _valueIndex = 1;

        /// <summary>
        /// ヘッダー行テキストを取得する。
        /// </summary>
        private static readonly string[] _headerLines = new[]
        {
            "",
            "########################################################################################################",
            "##",
            "## Lines beginning with # are comments.",
            "## Specify the ShortName or LongName of the option switch in Key and the setting value in Value.",
            "## If Value is left blank, it will be ignored.",
            "## Key and Value must be tab-separated.",
            "##",
            "########################################################################################################",
            "",
            "#Key\tValue"        };

        private readonly string _title;
        private readonly OptionCollection _options;

        /// <summary>
        /// パラメータファイルインスタンスを作成する。
        /// </summary>
        /// <param name="title">タイトル</param>
        /// <param name="options">オプション</param>
        public ParameterFile(string title, OptionCollection options)
        {
            _title = title;
            _options = options;
        }

        /// <summary>
        /// パラメータファイルからLongName-Value辞書を作成する。
        /// </summary>
        /// <param name="filePath">パラメータファイルPath</param>
        /// <returns>LongName-Value辞書</returns>
        public IReadOnlyDictionary<string, string> Parse(string filePath)
        {
            if (string.IsNullOrWhiteSpace(filePath)) throw new ArgumentNullException(nameof(filePath));
            if (!File.Exists(filePath)) throw new FileNotFoundException($"{filePath}が見つかりませんでした。");

            var paramsDict = new Dictionary<string, string>();

            using var reader = new StreamReader(filePath);
            string? line;
            while (!reader.EndOfStream)
            {
                line = reader.ReadLine();
                if (string.IsNullOrEmpty(line)) continue;
                if (line.StartsWith(_commentOut)) continue;

                var items = line.Split(_delimiter);
                if (items.Length == 1) continue;    // Keyのみはスキップ
                if (string.IsNullOrEmpty(items[_valueIndex])) continue;   // Valueがないのでスキップ

                var longName = _options.GetLongName(items[_keyIndex]);
                if (string.IsNullOrEmpty(longName)) throw new ArgumentException($"存在しないオプション名です。[{items[_keyIndex]}]");

                paramsDict[longName] = items[_valueIndex];
            }

            return paramsDict;
        }


        /// <summary>
        /// パラメータファイルを作成する。
        /// </summary>
        /// <param name="filePath">パラメータファイルPath</param>
        public void Create(string filePath)
        {
            if (string.IsNullOrWhiteSpace(filePath)) throw new ArgumentNullException(nameof(filePath));

            var lines = GetLines();
            using var writer = new StreamWriter(filePath);
            foreach (var line in lines)
            {
                writer.WriteLine(line);
            }
        }

        /// <summary>
        /// パラメータファイルに記載する行を取得する。
        /// </summary>
        /// <returns>パラメータファイルに記載する行排列</returns>
        private string[] GetLines()
        {
            return new[] { $"{_commentOut} {_title}" }
                .Concat(_headerLines)
                .Concat(_options.GetKeyValues())
                .ToArray();
        }
    }
}
